/* mathSelect.ts — select part of a typeset (StandardForm) result and copy it as input.
 *
 * The kernel's selectable StandardForm (print_latex.c, expr_to_latex_selectable) wraps every
 * subexpression in \htmlData{mp=K}{...}, which KaTeX turns into a span carrying data-mp="K",
 * and sends `sparts`, where sparts[K] is that subexpression's InputForm text. Typeset HTML has no
 * text worth copying -- it is positioned glyphs, plus a hidden MathML copy -- so the selection is
 * kept in terms of those spans instead, Mathematica-style:
 *
 *   click            selects the innermost subexpression under the pointer;
 *   click inside it  widens to the enclosing subexpression (repeat to climb to the whole result);
 *   drag             selects the smallest subexpression covering what was dragged over;
 *   Cmd/Ctrl+C       copies the selection's InputForm, ready to paste into an input cell;
 *   Escape           clears it.
 *
 * A widening step skips an enclosing span with the SAME text (a display wrapper such as HoldForm,
 * or an expression and its only rendered child), since selecting it would visibly do nothing; a
 * click once the whole result is selected starts again at the innermost subexpression.
 *
 * The highlight is a separate box over the selection, not a background on the span: the span is
 * inline, and an inline box is one line tall, so it would not cover a fraction or a radical. The box
 * is the union of the rectangles of everything the span contains -- except SVG: KaTeX draws a
 * radical's sign as a 400em-wide SVG clipped by its .hide-tail wrapper, and an element's rect
 * ignores the clip, so the wrapper (whose rect is the visible part) stands in for it.
 */

export interface MathSelectParams {
  /** sparts from the kernel; absent/empty → the action is inert (nothing selectable). */
  parts?: string[];
  /** Told the selected InputForm text (null when cleared), e.g. for a context-menu Copy. */
  onSelect?: (text: string | null) => void;
}

const NONE: string[] = [];

export function mathSelect(node: HTMLElement, params: MathSelectParams) {
  let parts: string[] = params.parts ?? NONE;
  let onSelect = params.onSelect;
  let sel: HTMLElement | null = null;
  let hl: HTMLDivElement | null = null;
  let dragged = false;          // the click that ends a drag must not also act as a click

  const active = () => parts.length > 0;

  function syncFocusable() {
    // Focusable so Cmd+C reaches THIS element (not an editor elsewhere) after a click.
    if (active()) node.tabIndex = -1; else node.removeAttribute('tabindex');
  }

  const textOf = (el: HTMLElement) => parts[Number(el.dataset.mp)] ?? '';

  /** The innermost tagged span at or above `el`, within this output. */
  function tagged(el: Element | null): HTMLElement | null {
    const t = el?.closest<HTMLElement>('[data-mp]') ?? null;
    return t && node.contains(t) ? t : null;
  }

  /** The nearest enclosing tagged span whose text differs from `el`'s. */
  function widen(el: HTMLElement): HTMLElement | null {
    const own = textOf(el);
    let p = tagged(el.parentElement);
    while (p && textOf(p) === own) p = tagged(p.parentElement);
    return p;
  }

  function drawHighlight() {
    if (!sel) { hl?.remove(); hl = null; return; }
    let l = Infinity, t = Infinity, r = -Infinity, b = -Infinity;
    for (const el of [sel, ...Array.from(sel.querySelectorAll('*'))]) {
      if (el.closest('svg')) continue;
      const rc = el.getBoundingClientRect();
      if (rc.width === 0 || rc.height === 0) continue;   // struts, spacing, empty vlists
      l = Math.min(l, rc.left); t = Math.min(t, rc.top);
      r = Math.max(r, rc.right); b = Math.max(b, rc.bottom);
    }
    if (!isFinite(l)) { hl?.remove(); hl = null; return; }
    const box = node.getBoundingClientRect();
    if (!hl) {
      hl = document.createElement('div');
      hl.className = 'mp-highlight';
      node.appendChild(hl);
    }
    /* Offsets inside `node` (position: relative), so the box scrolls with the output. */
    hl.style.left = `${l - box.left + node.scrollLeft - 1}px`;
    hl.style.top = `${t - box.top + node.scrollTop - 1}px`;
    hl.style.width = `${r - l + 2}px`;
    hl.style.height = `${b - t + 2}px`;
  }

  function setSel(el: HTMLElement | null) {
    sel = el;
    drawHighlight();
    onSelect?.(sel ? textOf(sel) : null);
  }

  function onClick(e: MouseEvent) {
    if (!active()) return;
    if (dragged) { dragged = false; return; }
    const hit = tagged(e.target as Element);
    if (!hit) { setSel(null); return; }
    setSel(sel && sel.isConnected && sel.contains(hit) ? (widen(sel) ?? hit) : hit);
  }

  function onMouseUp() {
    if (!active()) return;
    const ns = window.getSelection();
    if (!ns || ns.isCollapsed || ns.rangeCount === 0) return;
    const range = ns.getRangeAt(0);
    if (!node.contains(range.startContainer) || !node.contains(range.endContainer)) return;
    const common = range.commonAncestorContainer;
    const el = common.nodeType === Node.ELEMENT_NODE ? (common as Element) : common.parentElement;
    const hit = tagged(el);
    if (!hit) return;
    ns.removeAllRanges();        // show the snapped selection, not the raw glyph-run one
    dragged = true;
    setSel(hit);
  }

  function onCopy(e: ClipboardEvent) {
    if (!active() || !sel || !e.clipboardData) return;
    e.clipboardData.setData('text/plain', textOf(sel));
    e.preventDefault();
    const box = hl;                  // brief flash: the copy happened, and of what
    box?.classList.add('mp-flash');
    setTimeout(() => box?.classList.remove('mp-flash'), 450);
  }

  function onKeydown(e: KeyboardEvent) {
    if (e.key === 'Escape' && sel) { e.stopPropagation(); setSel(null); }
  }

  /* Reaching for another cell or output ends this selection, as in any editor. */
  function onFocusOut(e: FocusEvent) {
    if (!node.contains(e.relatedTarget as Node | null)) setSel(null);
  }

  syncFocusable();
  node.addEventListener('click', onClick);
  node.addEventListener('mouseup', onMouseUp);
  node.addEventListener('copy', onCopy);
  node.addEventListener('keydown', onKeydown);
  node.addEventListener('focusout', onFocusOut);

  return {
    update(p: MathSelectParams) {
      const next = p.parts ?? NONE;
      if (next !== parts) { parts = next; setSel(null); syncFocusable(); }
      onSelect = p.onSelect;
    },
    destroy() {
      node.removeEventListener('click', onClick);
      node.removeEventListener('mouseup', onMouseUp);
      node.removeEventListener('copy', onCopy);
      node.removeEventListener('keydown', onKeydown);
      node.removeEventListener('focusout', onFocusOut);
      hl?.remove();
    },
  };
}
