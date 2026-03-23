#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from html import escape
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from jinja2 import Environment, FileSystemLoader
import markdown as md_lib


DEFAULT_OUTPUT_DIR = "rendered_reports"
DEFAULT_CSS_PATH = "report_pdf.css"
MARKDOWN_EXTENSIONS = [
    "extra",
    "tables",
    "fenced_code",
    "toc",
    "sane_lists",
    "smarty",
]


class SampleIndexList(list):
    def __getitem__(self, item: Any) -> Any:
        if isinstance(item, str):
            for entry in self:
                if not isinstance(entry, dict):
                    continue
                sample_id = entry.get("collecting_lab_sample_id") or entry.get("sample_id")
                if sample_id == item:
                    return entry
            raise KeyError(item)
        return super().__getitem__(item)

    def get(self, item: str, default: Any = None) -> Any:
        try:
            return self[item]
        except KeyError:
            return default


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def is_lab_report_json(payload: Dict[str, Any]) -> bool:
    return isinstance(payload, dict) and "lab" in payload and "components" in payload


def get_lab_identifier(labdata: Dict[str, Any], fallback: str) -> str:
    lab_meta = labdata.get("lab", {})
    for key in ("lab_cod", "submitting_institution_id", "laboratory_name"):
        value = lab_meta.get(key)
        if value is not None and str(value).strip():
            return sanitize_filename(str(value).strip())
    return sanitize_filename(fallback)


def normalize_lab_payload(payload: Dict[str, Any], fallback_stem: str) -> Dict[str, Any]:
    if not is_lab_report_json(payload):
        return payload

    lab_meta = payload.setdefault("lab", {})
    if not lab_meta.get("lab_cod"):
        fallback_id = lab_meta.get("submitting_institution_id") or fallback_stem
        lab_meta["lab_cod"] = str(fallback_id)

    for component in payload.get("components", {}).values():
        if isinstance(component, dict) and "lab" not in component:
            component["lab"] = lab_meta
    return payload


def normalize_general_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    for component in payload.get("components", {}).values():
        if not isinstance(component, dict):
            continue
        for section_name in ("consensus", "variant", "typing", "qc", "metadata_metrics"):
            section = component.get(section_name)
            if isinstance(section, dict) and isinstance(section.get("samples"), list):
                section["samples"] = SampleIndexList(section["samples"])
    return payload


def sanitize_filename(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return safe.strip("._") or "report"


def collect_lab_jsons(labs_dir: Path) -> List[Tuple[Path, Dict[str, Any]]]:
    candidates: List[Tuple[Path, Dict[str, Any]]] = []
    for path in sorted(labs_dir.glob("*.json")):
        try:
            payload = load_json(path)
        except json.JSONDecodeError:
            continue
        if is_lab_report_json(payload):
            candidates.append((path, normalize_lab_payload(payload, path.stem)))
    return candidates


def safe_format_filter(value: Any, *args: Any, **kwargs: Any) -> str:
    if any(arg is None for arg in args) or any(v is None for v in kwargs.values()):
        return "N/A"
    try:
        return str(value) % (kwargs or args)
    except Exception:
        return "N/A"


def build_environment(template_path: Path) -> Environment:
    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        autoescape=False,
        keep_trailing_newline=True,
    )
    env.filters["format"] = safe_format_filter
    return env


def render_template(
    template_path: Path,
    general_data: Dict[str, Any],
    labdata: Optional[Dict[str, Any]],
) -> str:
    env = build_environment(template_path)
    template = env.get_template(template_path.name)
    rendered = template.render(general=general_data, labdata=labdata)
    return postprocess_rendered_markdown(rendered)


def normalize_table_spacing(markdown_text: str) -> str:
    lines = markdown_text.splitlines()
    normalized: List[str] = []
    total = len(lines)

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("|"):
            normalized.append(stripped)
            continue

        if stripped == "":
            prev_is_table = bool(normalized and normalized[-1].strip().startswith("|"))
            next_line = lines[idx + 1].strip() if idx + 1 < total else ""
            next_is_table = next_line.startswith("|")
            if prev_is_table and next_is_table:
                continue

        normalized.append(line)

    return "\n".join(normalized)


def normalize_figure_blocks(markdown_text: str) -> str:
    markdown_text = re.sub(r"\n[ \t]+(<(?:figure|img|figcaption|/figure)[^>]*>)", r"\n\1", markdown_text)
    markdown_text = re.sub(r"(?m)^[ \t]*(<figure>)$", r"\1", markdown_text)
    markdown_text = re.sub(r"(?m)^[ \t]*(</figure>)$", r"\1", markdown_text)
    markdown_text = re.sub(r"(?m)^[ \t]*(<img\b[^>]*>)$", r"\1", markdown_text)
    markdown_text = re.sub(r"(?m)^[ \t]*(<figcaption>.*?</figcaption>)$", r"\1", markdown_text)
    return markdown_text


def latex_block_to_html(expression: str) -> str:
    compact = " ".join(expression.split())
    fraction_match = re.match(
        r"\\text\{(?P<lhs>.+?)\}\s*=\s*\\frac\{\s*\\text\{(?P<num>.+?)\}\s*\}\s*\{\s*\\text\{(?P<den>.+?)\}\s*\}",
        compact,
    )
    if fraction_match:
        lhs = escape(fraction_match.group("lhs"))
        numerator = escape(fraction_match.group("num"))
        denominator = escape(fraction_match.group("den"))
        return (
            '<div class="equation">'
            f'<span class="equation-lhs">{lhs}</span>'
            '<span class="equation-equals">=</span>'
            '<span class="equation-fraction">'
            f'<span class="equation-numerator">{numerator}</span>'
            f'<span class="equation-denominator">{denominator}</span>'
            "</span>"
            "</div>"
        )

    fallback = escape(expression.strip())
    return f'<div class="equation"><code>{fallback}</code></div>'


def replace_display_math_blocks(markdown_text: str) -> str:
    return re.sub(
        r"\$\$\s*(.*?)\s*\$\$",
        lambda match: "\n\n" + latex_block_to_html(match.group(1)) + "\n\n",
        markdown_text,
        flags=re.DOTALL,
    )


def normalize_missing_markers(markdown_text: str) -> str:
    markdown_text = re.sub(r"(?<![\w/.-])None(?![\w/.-])", "N/A", markdown_text)
    markdown_text = re.sub(r"(?m)(\|\s*)NA(\s*(?=\|))", r"\1N/A\2", markdown_text)
    return markdown_text


def postprocess_rendered_markdown(markdown_text: str) -> str:
    cleaned = markdown_text.lstrip()
    cleaned = normalize_figure_blocks(cleaned)
    cleaned = normalize_table_spacing(cleaned)
    cleaned = replace_display_math_blocks(cleaned)
    cleaned = normalize_missing_markers(cleaned)
    return cleaned


def postprocess_rendered_html(html_text: str) -> str:
    html_text = re.sub(r"<pre><code>\s*(<figure>.*?</figure>)\s*</code></pre>", r"\1", html_text, flags=re.DOTALL)
    html_text = re.sub(r"<p>\s*(<figure>.*?</figure>)\s*</p>", r"\1", html_text, flags=re.DOTALL)
    html_text = re.sub(r"<p>\s*(<div class=\"equation\">.*?</div>)\s*</p>", r"\1", html_text, flags=re.DOTALL)
    return html_text


def markdown_to_html(markdown_text: str, title: str, css_text: str, base_dir: Path) -> str:
    html_body = md_lib.markdown(markdown_text, extensions=MARKDOWN_EXTENSIONS)
    html_body = postprocess_rendered_html(html_body)
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{escape(title)}</title>
    <base href="{base_dir.resolve().as_uri()}/">
    <style>
{css_text}
    </style>
  </head>
  <body>
    {html_body}
  </body>
</html>
"""


def detect_browser_backend() -> Optional[str]:
    for candidate in ("google-chrome", "chromium", "chromium-browser"):
        executable = shutil.which(candidate)
        if executable:
            return executable
    return None


def render_pdf_with_chrome(html_path: Path, pdf_path: Path, chrome_path: str) -> None:
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        chrome_path,
        "--headless=new",
        "--disable-gpu",
        "--no-pdf-header-footer",
        "--allow-file-access-from-files",
        f"--print-to-pdf={pdf_path.resolve()}",
        html_path.resolve().as_uri(),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "PDF generation failed with Chrome.\n"
            f"Command: {' '.join(command)}\n"
            f"stderr: {result.stderr.strip()}"
        )


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def discover_existing_markdown_reports(markdown_root: Path) -> List[Dict[str, Any]]:
    reports: List[Dict[str, Any]] = []
    for path in sorted(markdown_root.rglob("*.md")):
        rel_parent = path.relative_to(markdown_root).parent
        stem = path.stem
        reports.append(
            {
                "kind": "existing_markdown",
                "identifier": stem,
                "title": stem.replace("_", " ").strip() or stem,
                "markdown_text": path.read_text(encoding="utf-8"),
                "subdir": rel_parent,
                "stem": stem,
                "source_markdown_path": path,
            }
        )
    return reports


def build_report_targets(
    general_data: Dict[str, Any],
    template_path: Path,
    labs_dir: Optional[Path],
) -> List[Dict[str, Any]]:
    reports: List[Dict[str, Any]] = []
    reports.append(
        {
            "kind": "general",
            "identifier": "general",
            "title": "RELECOV 2026 Dry-Lab EQA General Report",
            "markdown_text": render_template(template_path, general_data, labdata=None),
            "subdir": Path(),
            "stem": "general_report",
        }
    )

    if not labs_dir:
        return reports

    lab_entries = collect_lab_jsons(labs_dir)
    if not lab_entries:
        raise FileNotFoundError(
            f"No compatible individual lab JSON files were found in {labs_dir}. "
            "Expected files with top-level 'lab' and 'components' keys."
        )

    for path, payload in lab_entries:
        lab_id = get_lab_identifier(payload, path.stem)
        reports.append(
            {
                "kind": "lab",
                "identifier": lab_id,
                "title": f"RELECOV 2026 Dry-Lab EQA Technical Report - {lab_id}",
                "markdown_text": render_template(template_path, general_data, labdata=payload),
                "subdir": Path("labs"),
                "stem": f"lab_{lab_id}",
            }
        )

    return reports


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render the RELECOV report Jinja markdown template into a general report and, "
            "optionally, individual laboratory reports. The rendered markdown can also be "
            "exported to PDF using a custom CSS layout."
        )
    )
    parser.add_argument("--template", required=False, help="Path to the markdown Jinja template.")
    parser.add_argument("--general-json", required=False, help="Path to general.json.")
    parser.add_argument(
        "--labs-dir",
        default=None,
        help=(
            "Optional directory containing individual lab JSON files compatible with the template "
            "(for example lab_*.json). If provided, individual lab reports are rendered too."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory where rendered reports will be written. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--css",
        default=DEFAULT_CSS_PATH,
        help=f"CSS file used for PDF export. Default: {DEFAULT_CSS_PATH}",
    )
    parser.add_argument(
        "--markdown-only",
        action="store_true",
        help="Render markdown outputs only and skip PDF export.",
    )
    parser.add_argument(
        "--pdf-only",
        action="store_true",
        help="Export PDFs only. Markdown is rendered in memory but not saved.",
    )
    parser.add_argument(
        "--keep-html",
        action="store_true",
        help="Keep the intermediate HTML files used for PDF generation.",
    )
    parser.add_argument(
        "--input-markdown-dir",
        default=None,
        help=(
            "Optional directory containing already-rendered markdown reports. "
            "If provided together with PDF export, PDFs are generated directly from those markdown files "
            "without re-rendering the Jinja template."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    if args.markdown_only and args.pdf_only:
        raise SystemExit("Use either --markdown-only or --pdf-only, not both.")

    render_markdown = not args.pdf_only
    render_pdf = not args.markdown_only

    template_path = Path(args.template).resolve() if args.template else None
    general_json_path = Path(args.general_json).resolve() if args.general_json else None
    labs_dir = Path(args.labs_dir).resolve() if args.labs_dir else None
    input_markdown_dir = Path(args.input_markdown_dir).resolve() if args.input_markdown_dir else None
    output_dir = Path(args.output_dir).resolve()
    css_path = Path(args.css).resolve()

    render_from_existing_markdown = input_markdown_dir is not None

    if render_from_existing_markdown and render_markdown:
        raise SystemExit("--input-markdown-dir can only be used with PDF export. Use --pdf-only or omit --markdown-only.")

    if not render_from_existing_markdown:
        if template_path is None or not template_path.exists():
            raise SystemExit(f"Template not found: {template_path}")
        if general_json_path is None or not general_json_path.exists():
            raise SystemExit(f"General JSON not found: {general_json_path}")
    if labs_dir and not labs_dir.exists():
        raise SystemExit(f"Labs directory not found: {labs_dir}")
    if input_markdown_dir and not input_markdown_dir.exists():
        raise SystemExit(f"Input markdown directory not found: {input_markdown_dir}")
    if not css_path.exists():
        raise SystemExit(f"CSS file not found: {css_path}")

    css_text = css_path.read_text(encoding="utf-8")
    if render_from_existing_markdown:
        reports = discover_existing_markdown_reports(input_markdown_dir)
        if not reports:
            raise SystemExit(f"No markdown files were found under {input_markdown_dir}")
        base_dir = input_markdown_dir
    else:
        general_data = normalize_general_payload(load_json(general_json_path))
        reports = build_report_targets(general_data, template_path, labs_dir)
        base_dir = template_path.parent

    markdown_root = output_dir / "markdown"
    pdf_root = output_dir / "pdf"
    html_root = output_dir / "html"

    browser_backend = None
    if render_pdf:
        browser_backend = detect_browser_backend()
        if browser_backend is None:
            raise SystemExit(
                "PDF export requires a Chrome/Chromium executable available in PATH. "
                "None of google-chrome, chromium, or chromium-browser was found."
            )

    rendered_markdown_paths: List[Tuple[Dict[str, Any], Path]] = []

    for report in reports:
        markdown_path = markdown_root / report["subdir"] / f"{report['stem']}.md"
        rendered_markdown_paths.append((report, markdown_path))
        if render_markdown:
            write_text(markdown_path, report["markdown_text"])

    if render_pdf:
        for report, markdown_path in rendered_markdown_paths:
            html_path = html_root / report["subdir"] / f"{report['stem']}.html"
            pdf_path = pdf_root / report["subdir"] / f"{report['stem']}.pdf"
            html_text = markdown_to_html(report["markdown_text"], report["title"], css_text, base_dir)
            write_text(html_path, html_text)
            render_pdf_with_chrome(html_path, pdf_path, browser_backend)
            if not args.keep_html:
                html_path.unlink(missing_ok=True)

    print(f"Rendered {len(reports)} report(s).")
    if render_markdown:
        print(f"Markdown output: {markdown_root}")
    if render_pdf:
        print(f"PDF output: {pdf_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
