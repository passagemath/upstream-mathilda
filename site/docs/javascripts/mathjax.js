window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  /* Emit the CSS for every font glyph up front, as static text in MathJax's
   * <style id="MJX-CHTML-styles">, instead of adding one rule per glyph on
   * demand via CSSOM insertRule (the default "adaptive" mode). Material's
   * instant navigation (`navigation.instant`) swaps pages without a reload and
   * can leave that <style> detached; adaptive rules then fail to insert (its
   * `sheet` is null) and every glyph the first page did not use renders blank
   * -- most of the formulae on any page reached by a link. Static text
   * survives detaching and re-attaching; CSSOM-inserted rules do not. */
  chtml: {
    adaptiveCSS: false
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

/* Re-typeset after every page load, including instant navigation. */
document$.subscribe(() => {
  if (typeof MathJax === "undefined" || !MathJax.typesetPromise)
    return;   /* first load: MathJax's own startup typesets the page */
  const styles = MathJax.startup.output && MathJax.startup.output.chtmlStyles;
  if (styles && !styles.isConnected) document.head.appendChild(styles);
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
