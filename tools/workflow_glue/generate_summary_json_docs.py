"""Render the public summary.json reference embedded in README.md."""

import argparse
import json
from pathlib import Path


START = "<!-- BEGIN GENERATED SUMMARY.JSON REFERENCE -->"
END = "<!-- END GENERATED SUMMARY.JSON REFERENCE -->"


HIDDEN_TOP_LEVEL_KEYS = frozenset(("schema_id", "schema_version"))


def argparser():
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--schema", required=True, type=Path, help="summary.json schema")
    parser.add_argument("--readme", required=True, type=Path, help="README containing generated markers")
    parser.add_argument("--check", action="store_true", help="fail if the README block is stale")
    return parser


def ref_name(schema):
    """Return a compact display name for a local schema reference."""
    if "allOf" in schema:
        for part in schema["allOf"]:
            if "$ref" in part:
                return ref_name(part)
    reference = schema.get("$ref")
    if reference and reference.startswith("#/$defs/"):
        return f"`{reference.rsplit('/', 1)[1]}`"
    if "const" in schema:
        return f"constant `{schema['const']}`"
    if "enum" in schema:
        return "one of " + ", ".join(f"`{value}`" for value in schema["enum"])
    data_type = schema.get("type")
    if isinstance(data_type, list):
        return " or ".join(f"`{value}`" for value in data_type)
    if data_type:
        return f"`{data_type}`"
    if "anyOf" in schema:
        return "optional value (`null` when unavailable)"
    if "allOf" in schema:
        return "composed object"
    return "value"


def description(schema):
    """Return required public documentation prose."""
    try:
        return schema["description"]
    except KeyError as error:
        raise ValueError("A documented summary.json field is missing a description") from error


def title(schema):
    """Return a required human-facing schema title."""
    try:
        return schema["title"]
    except KeyError as error:
        raise ValueError("A documented top-level summary.json field is missing a title") from error


def table(rows):
    lines = ["| Field | Type | Description |", "| --- | --- | --- |"]
    lines.extend(f"| `{name}` | {kind} | {text} |" for name, kind, text in rows)
    return "\n".join(lines)


def object_rows(schema):
    return [
        (name, ref_name(value), description(value))
        for name, value in schema.get("properties", {}).items()
    ]


def column_rows(schema):
    rows = []
    for column in schema.get("prefixItems", []):
        name = column.get("properties", {}).get("name", {}).get("const")
        if name is None:
            continue
        name_schema = column["properties"]["name"]
        rows.append((name, ref_name(column), description(name_schema)))
    if schema.get("items") not in (None, False):
        pattern = schema["items"].get("properties", {}).get("name", {}).get("pattern")
        if pattern:
            rows.append((pattern, ref_name(schema["items"]), description(schema["items"])))
    return rows


def resolve(schema, defs):
    """Follow local references and compose the relevant non-null schema branch."""
    while "$ref" in schema:
        reference = schema["$ref"]
        if not reference.startswith("#/$defs/"):
            raise ValueError(f"Unsupported schema reference: {reference}")
        schema = defs[reference.rsplit("/", 1)[1]]
    if "anyOf" in schema:
        for option in schema["anyOf"]:
            if option.get("type") != "null":
                schema = resolve(option, defs)
                break
        else:
            raise ValueError("Optional summary.json value has no non-null schema")

    if "allOf" in schema:
        combined = {}
        for part in schema["allOf"]:
            part = resolve(part, defs)
            for key, value in part.items():
                if key == "properties":
                    combined.setdefault(key, {}).update(value)
                elif key == "required":
                    combined.setdefault(key, []).extend(value)
                elif key != "allOf":
                    combined[key] = value
        for key, value in schema.items():
            if key != "allOf" and key not in combined:
                combined[key] = value
        return combined
    return schema


def dataframe_column_layouts(schema, defs, path=""):
    """Yield dataframe paths and their resolved column schemas within a section."""
    schema = resolve(schema, defs)
    properties = schema.get("properties", {})
    kind = properties.get("kind", {})
    if kind.get("const") == "dataframe":
        columns = properties.get("columns")
        if columns is None:
            raise ValueError("Dataframe schema is missing columns")
        yield path, resolve(columns, defs)
        return

    if schema.get("type") == "array" and "items" in schema:
        yield from dataframe_column_layouts(schema["items"], defs, f"{path}[]")
        return

    for name, child in properties.items():
        child_path = f"{path}.{name}" if path else name
        yield from dataframe_column_layouts(child, defs, child_path)


def layout_label(root_name, paths):
    """Return concise inline Markdown for one or more dataframe paths."""
    display_paths = [path or root_name for path in paths]
    quoted = [f"`{path}`" for path in display_paths]
    if len(quoted) == 1:
        return quoted[0]
    if len(quoted) == 2:
        return " and ".join(quoted)
    return ", ".join(quoted[:-1]) + f", and {quoted[-1]}"


def grouped_dataframe_layouts(schema, defs):
    """Group dataframe paths that use exactly the same column layout."""
    groups = {}
    for path, columns in dataframe_column_layouts(schema, defs):
        key = json.dumps(columns, sort_keys=True)
        if key not in groups:
            groups[key] = (columns, [])
        groups[key][1].append(path)
    return groups.values()


def root_definition(schema, defs):
    """Return the object schema documented for a top-level property."""
    schema = resolve(schema, defs)
    if schema.get("type") == "array" and "items" in schema:
        return resolve(schema["items"], defs)
    return schema


def render(schema):
    """Render the generated per-top-level-key content reference."""
    defs = schema["$defs"]
    lines = [START]

    for root_name, root in schema["properties"].items():
        if root_name in HIDDEN_TOP_LEVEL_KEYS:
            continue
        definition = root_definition(root, defs)
        lines.extend([f"\n#### `{root_name}` — {title(root)}\n", description(root), ""])
        rows = object_rows(definition)
        if rows:
            lines.append(table(rows))
        else:
            lines.append(f"**Type:** {ref_name(root)}.")
        for columns, paths in grouped_dataframe_layouts(definition, defs):
            lines.extend([
                f"\nThe **{layout_label(root_name, paths)}** dataframe columns are:\n",
                table(column_rows(columns)),
                "",
            ])

    lines.append(END)
    return "\n".join(lines)


def replace_block(readme, generated):
    try:
        before, rest = readme.split(START, 1)
        _, after = rest.split(END, 1)
    except ValueError as error:
        raise ValueError("README is missing summary.json reference markers") from error
    return before + generated + after


def main(args):
    schema = json.loads(args.schema.read_text(encoding="utf-8"))
    generated = render(schema)
    current = args.readme.read_text(encoding="utf-8")
    updated = replace_block(current, generated)
    if args.check:
        if updated != current:
            raise SystemExit(f"{args.readme} summary.json reference is out of date")
        return 0
    args.readme.write_text(updated, encoding="utf-8")
    return 0
