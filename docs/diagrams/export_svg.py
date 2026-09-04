"""Export standalone light/dark SVGs from the delivered archify HTML.

The interactive HTML (docs/workflow.html) styles its inline SVG through CSS
custom properties that switch with [data-theme]. GitHub renders SVG inside
<img>, where scripts are blocked and some rasterizers ignore var(), so this
script copies the SVG out, keeps only the rules that target SVG classes, and
substitutes each theme's concrete variable values before writing:

    assets/eca-pp-workflow-light.svg
    assets/eca-pp-workflow-dark.svg

Usage (from the repository root):
    python docs/diagrams/export_svg.py
"""

from __future__ import annotations

import re
import sys
import xml.dom.minidom
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HTML = ROOT / "docs" / "workflow.html"
ASSETS = ROOT / "assets"
THEMES = ("light", "dark")


def split_rules(css: str) -> list[tuple[str, str]]:
    """Top-level (selector, body) pairs; at-rules (@media, @keyframes) are dropped."""
    rules = []
    i, n = 0, len(css)
    while i < n:
        j = css.find("{", i)
        if j < 0:
            break
        selector = re.sub(r"/\*.*?\*/", "", css[i:j], flags=re.DOTALL).strip()
        depth, k = 1, j + 1
        while k < n and depth:
            depth += {"{": 1, "}": -1}.get(css[k], 0)
            k += 1
        if not selector.startswith("@"):
            rules.append((selector, css[j + 1 : k - 1]))
        i = k
    return rules


def variables(body: str) -> dict[str, str]:
    return {m.group(1): m.group(2).strip() for m in re.finditer(r"(--[\w-]+)\s*:\s*([^;]+);", body)}


def theme_variables(rules: list[tuple[str, str]], theme: str) -> dict[str, str]:
    base, override = {}, {}
    for selector, body in rules:
        parts = [p.strip() for p in selector.split(",")]
        if ":root" in parts:
            base.update(variables(body))
        if parts == [f'[data-theme="{theme}"]']:
            override.update(variables(body))
    resolved = {**base, **override}
    for _ in range(8):  # variables may reference other variables
        changed = False
        for key, value in resolved.items():
            new = re.sub(r"var\((--[\w-]+)(?:,\s*([^)]+))?\)",
                         lambda m: resolved.get(m.group(1), m.group(2) or ""), value)
            if new != value:
                resolved[key], changed = new, True
        if not changed:
            break
    return resolved


def substitute(css: str, values: dict[str, str]) -> str:
    return re.sub(r"var\((--[\w-]+)(?:,\s*([^)]+))?\)",
                  lambda m: values.get(m.group(1), m.group(2) or "inherit"), css)


def keep_rule(selector: str, svg_classes: set[str]) -> bool:
    if "data-theme" in selector or ":root" in selector or "#theme-icon" in selector:
        return False
    if "data-preset" in selector and 'data-preset="classic"' not in selector:
        return False
    for part in (p.strip() for p in selector.split(",")):
        if "html" in part or "body" in part:
            continue
        if any(f".{c}" in part for c in svg_classes):
            return True
        if re.search(r"(^|[\s>])(svg|text|marker|path|rect)(\b|[.:\[])", part):
            return True
    return False


def main() -> int:
    html = HTML.read_text(encoding="utf-8")
    style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)
    svg_match = re.search(r"<svg.*?</svg>", html, re.DOTALL)
    if style_match is None or svg_match is None:
        print(f"error: no <style> or <svg> block found in {HTML}", file=sys.stderr)
        return 1
    style, svg = style_match.group(1), svg_match.group(0)
    rules = split_rules(style)
    svg_classes = {c for m in re.findall(r'class="([^"]*)"', svg) for c in m.split()}
    kept = [(s, b) for s, b in rules if keep_rule(s, svg_classes)]
    missing = [c for c in sorted(svg_classes) if not any(f".{c}" in s for s, _ in kept)]
    if missing:
        print(f"warning: no CSS rule kept for classes {missing}", file=sys.stderr)

    head = svg[: svg.index(">") + 1]  # the opening <svg ...> tag
    body = svg[len(head):]
    ASSETS.mkdir(exist_ok=True)
    for theme in THEMES:
        values = theme_variables(rules, theme)
        compact = ((s, re.sub(r"\s+", " ", b).strip()) for s, b in kept)
        css = "\n".join(f"{s}{{{b}}}" for s, b in compact)
        css = substitute(css, values)
        css = re.sub(r"/\*.*?\*/", "", css, flags=re.DOTALL)
        background = values.get("--bg", "#ffffff")
        new_head = head.replace(
            "<svg ", f'<svg xmlns="http://www.w3.org/2000/svg" data-theme="{theme}" ', 1
        )
        # A real rect, not a root style attribute: rasterizers such as cairosvg
        # ignore CSS backgrounds on the outermost <svg>.
        backdrop = f'<rect width="100%" height="100%" fill="{background}"/>'
        out = ASSETS / f"eca-pp-workflow-{theme}.svg"
        out.write_text(
            f"{new_head}\n<style><![CDATA[\n{css}\n]]></style>\n{backdrop}{body}", encoding="utf-8"
        )
        xml.dom.minidom.parse(str(out))  # raises if not well-formed
        print(f"{out.relative_to(ROOT)}: {out.stat().st_size} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
