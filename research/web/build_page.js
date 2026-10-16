/* Assemble dist/ladder.html: inline the Observatory solver, the page data and the app, and
 * pre-render the equations to SVG with MathJax (no runtime math library).
 * Usage: node build_page.js   (run build_data.js first) */
"use strict";
const fs = require("fs"), path = require("path");
const W = path.join(__dirname, "../../web");
const req = (m) => require(path.join(W, "node_modules", m));
req("mathjax-full/js/util/asyncLoad/node.js");
const { mathjax } = req("mathjax-full/js/mathjax.js");
const { TeX } = req("mathjax-full/js/input/tex.js");
const { SVG } = req("mathjax-full/js/output/svg.js");
const { liteAdaptor } = req("mathjax-full/js/adaptors/liteAdaptor.js");
const { RegisterHTMLHandler } = req("mathjax-full/js/handlers/html.js");
const { AllPackages } = req("mathjax-full/js/input/tex/AllPackages.js");
const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const doc = mathjax.document("", { InputJax: new TeX({ packages: AllPackages }), OutputJax: new SVG({ fontCache: "none" }) });
const r = String.raw;
const EQ = {
  word: r`w_k \;=\; (\mathtt{Ba})^{k+1}\,(\mathtt{bA})^{k+1}\;\;\sim\;\;\mathtt{A}\,(\mathtt{Ba})^{k}\,\mathtt{B}\,\mathtt{a}\,(\mathtt{bA})^{k}\,\mathtt{b},\qquad \lvert w_k\rvert = 4k+4`,
  slice: r`\mathbf r_1=(-1,0),\;\; \mathbf r_2=(1,0),\;\; \mathbf r_3=(0,0),\qquad \dot{\mathbf r}_1=\dot{\mathbf r}_2=(v_1,v_2),\;\; \dot{\mathbf r}_3=-\frac{2}{m_3}\,(v_1,v_2),\qquad T^{*}=T\,\lvert E\rvert^{3/2}`,
};
function tex(src) {
  const html = adaptor.outerHTML(doc.convert(src, { display: true }));
  return `<div class="eq">${html.replace("<mjx-container", `<mjx-container role="math" aria-label="${src.replace(/\s+/g, " ").replace(/"/g, "&quot;")}"`)}</div>`;
}
let page = fs.readFileSync(path.join(__dirname, "src/ladder.html"), "utf8");
page = page.replace(/\{\{EQ:(\w+)\}\}/g, (_, k) => { if (!EQ[k]) throw new Error("missing EQ " + k); return tex(EQ[k]); });
const mjx = `mjx-container { display: inline-block; } mjx-container[display="true"] { display: block; margin: 0; } mjx-container svg { overflow: visible; }`;
page = page.split("</style>").join(mjx + "\n</style>");
const data = fs.readFileSync(path.join(__dirname, "ladder_data.json"), "utf8");
page = page.split("/*__SOLVER__*/").join(fs.readFileSync(path.join(W, "solver.js"), "utf8"));
page = page.split("/*__DATA__*/null").join(data);
page = page.split("/*__APP__*/").join(fs.readFileSync(path.join(__dirname, "src/ladder.js"), "utf8"));
fs.mkdirSync(path.join(__dirname, "dist"), { recursive: true });
const out = path.join(__dirname, "dist/ladder.html");
fs.writeFileSync(out, page);
console.log(out, (page.length / 1024).toFixed(0) + " KB");
// ladder.html is a body fragment for the artifact host, which adds the document head; this copy opens
// directly in a browser
const standalone = path.join(__dirname, "dist/ladder.standalone.html");
fs.writeFileSync(standalone, `<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style>\n</head>\n<body>\n${page}\n</body>\n</html>\n`);
console.log(standalone);
