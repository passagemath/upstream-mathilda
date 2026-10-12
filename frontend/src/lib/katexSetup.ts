/* katexSetup.ts — teach KaTeX the StandardForm characters it lacks: ⅇ ⅈ ⅆ and ≠.
 *
 * The kernel's StandardForm LaTeX (print_latex.c) spells the constants E and I
 * as ⅇ (U+2147) and ⅈ (U+2148), and an integral's differential as ⅆ (U+2146),
 * exactly as Mathematica's StandardForm does -- it is how an upright-everything
 * notation keeps the constant E apart from a variable called e.
 *
 * KaTeX ships no glyphs or metrics for the Letterlike Symbols block, so out of
 * the box each one costs a console warning per render ("Unrecognized Unicode
 * character", then "No character metrics") and is laid out as a zero-height box,
 * which mis-seats a superscript on it (ⅇ^x). Both are fixed through KaTeX's own
 * extension hooks rather than by patching the library:
 *
 *   - __defineSymbol makes each a known math-mode symbol (an ordinary, set in
 *     the Main font), so the parser no longer treats it as foreign text;
 *   - __setFontMetrics gives Main-Regular an entry for each, borrowed from the
 *     plain letter it doubles (d, e, i). The glyph itself comes from the
 *     browser's font fallback; only its box comes from these numbers, and the
 *     plain letter's box is the right size for it.
 *
 * Import for its side effect, once, before the first render. Idempotent.
 */
import katex from 'katex';
import fontMetricsData from 'katex/src/fontMetricsData.js';

type Metrics = Record<string, Record<string, number[]>>;

/** [double-struck code point, the plain letter whose metrics it borrows] */
const DOUBLE_STRUCK: [number, string][] = [
  [0x2146, 'd'],   // ⅆ  differential d
  [0x2147, 'e'],   // ⅇ  exponential e (the constant E)
  [0x2148, 'i'],   // ⅈ  imaginary i (the constant I)
];

/* ≠ as ONE glyph. StandardForm writes Unequal as \char"2260 rather than \neq,
   because KaTeX composes \neq from `=` and a slash out of its own font, and that
   slash does not exist in the Arial StandardForm is set in. \char takes the
   character as-is, so it only needs a box: the `=` it negates has the right one.
   (Not a defined symbol: KaTeX already maps a typed ≠ to its composite \neq.) */
const BORROWED_ONLY: [number, string][] = [
  [0x2260, '='],   // ≠
];

/* The hooks are KaTeX's extension surface (mhchem uses them) but are not in its
   published typings, hence the narrow cast. */
const k = katex as unknown as {
  __defineSymbol(mode: string, font: string, group: string, replace: string, name: string): void;
  __setFontMetrics(font: string, metrics: Record<string, number[]>): void;
};

const main = { ...(fontMetricsData as Metrics)['Main-Regular'] };
for (const [cp, letter] of DOUBLE_STRUCK) {
  const ch = String.fromCodePoint(cp);
  main[String(cp)] = main[String(letter.charCodeAt(0))];
  k.__defineSymbol('math', 'main', 'textord', ch, ch);
}
for (const [cp, ch] of BORROWED_ONLY) main[String(cp)] = main[String(ch.charCodeAt(0))];
k.__setFontMetrics('Main-Regular', main);
