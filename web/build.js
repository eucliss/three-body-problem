/* Assemble dist/index.html: inline the solver, the app, the refined periodic orbits and the
 * fate map, and pre-render every equation to SVG with MathJax (no runtime math library). */
const fs = require("fs");
const path = require("path");
require("mathjax-full/js/util/asyncLoad/node.js");
const { mathjax } = require("mathjax-full/js/mathjax.js");
const { TeX } = require("mathjax-full/js/input/tex.js");
const { SVG } = require("mathjax-full/js/output/svg.js");
const { liteAdaptor } = require("mathjax-full/js/adaptors/liteAdaptor.js");
const { RegisterHTMLHandler } = require("mathjax-full/js/handlers/html.js");
const { AllPackages } = require("mathjax-full/js/input/tex/AllPackages.js");

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const doc = mathjax.document("", { InputJax: new TeX({ packages: AllPackages }), OutputJax: new SVG({ fontCache: "none" }) });
const r = String.raw;

const EQ = {
  hero: r`\ddot{\mathbf r}_i=-G\sum_{j\neq i}\frac{m_j\,(\mathbf r_i-\mathbf r_j)}{\lvert\mathbf r_i-\mathbf r_j\rvert^{3}}`,
  main: r`\ddot{\mathbf r}_i \;=\; -\,G\sum_{j\neq i} m_j\,\frac{\mathbf r_i-\mathbf r_j}{\left(\lvert\mathbf r_i-\mathbf r_j\rvert^{2}+\varepsilon^{2}\right)^{3/2}},\qquad i=1,2,3`,
  expanded: r`\begin{aligned}
\ddot{\mathbf r}_1 &= -\,G m_2\,\frac{\mathbf r_1-\mathbf r_2}{r_{12}^{3}} \;-\; G m_3\,\frac{\mathbf r_1-\mathbf r_3}{r_{13}^{3}}\\[4pt]
\ddot{\mathbf r}_2 &= -\,G m_3\,\frac{\mathbf r_2-\mathbf r_3}{r_{23}^{3}} \;-\; G m_1\,\frac{\mathbf r_2-\mathbf r_1}{r_{12}^{3}}\\[4pt]
\ddot{\mathbf r}_3 &= -\,G m_1\,\frac{\mathbf r_3-\mathbf r_1}{r_{13}^{3}} \;-\; G m_2\,\frac{\mathbf r_3-\mathbf r_2}{r_{23}^{3}}
\end{aligned}\qquad r_{ij}=\lvert\mathbf r_i-\mathbf r_j\rvert`,
  hamiltonian: r`H=\sum_{i=1}^{3}\frac{\lvert\mathbf p_i\rvert^{2}}{2m_i}-\sum_{i<j}\frac{G\,m_i m_j}{r_{ij}},\qquad
\dot{\mathbf r}_i=\frac{\partial H}{\partial\mathbf p_i},\quad \dot{\mathbf p}_i=-\frac{\partial H}{\partial\mathbf r_i}`,
  symmetric: r`\ddot{\mathbf s}_i=-\,G M\,\frac{\mathbf s_i}{s_i^{3}}+G\,m_i\,\mathbf S,\qquad
\mathbf S=\frac{\mathbf s_1}{s_1^{3}}+\frac{\mathbf s_2}{s_2^{3}}+\frac{\mathbf s_3}{s_3^{3}},\qquad M=m_1+m_2+m_3`,
  scaling: r`t_*=\sqrt{\frac{L^{3}}{G M}},\qquad \mathbf r=L\,\tilde{\mathbf r},\;\; t=t_*\,\tilde t,\;\; m_i=M\,\tilde m_i
\;\;\Longrightarrow\;\; \ddot{\tilde{\mathbf r}}_i=-\sum_{j\neq i}\tilde m_j\,\frac{\tilde{\mathbf r}_i-\tilde{\mathbf r}_j}{\lvert\tilde{\mathbf r}_i-\tilde{\mathbf r}_j\rvert^{3}}`,
  integrals: r`\begin{aligned}
E&=\sum_i \tfrac12 m_i\lvert\dot{\mathbf r}_i\rvert^{2}-\sum_{i<j}\frac{G m_i m_j}{r_{ij}} &&\text{energy}\\
\mathbf P&=\sum_i m_i\,\dot{\mathbf r}_i &&\text{momentum}\\
\mathbf L&=\sum_i m_i\,\mathbf r_i\times\dot{\mathbf r}_i &&\text{angular momentum}\\
\mathbf C&=\sum_i m_i\,\mathbf r_i-\mathbf P\,t &&\text{centre-of-mass drift}
\end{aligned}`,
  sundman: r`\mathbf r_i(t)=\sum_{k=0}^{\infty}\mathbf a_{i,k}\,\sigma^{k},\qquad \lvert\sigma\rvert<1\ \text{ covers }\ -\infty<t<\infty`,
  lagomega: r`\omega^{2}=\frac{G\,(m_1+m_2+m_3)}{a^{3}}`,
  routh: r`27\,(m_1m_2+m_2m_3+m_3m_1)\;<\;(m_1+m_2+m_3)^{2}`,
  quintic: r`(m_1+m_2)\lambda^{5}+(3m_1+2m_2)\lambda^{4}+(3m_1+m_2)\lambda^{3}-(m_2+3m_3)\lambda^{2}-(2m_2+3m_3)\lambda-(m_2+m_3)=0`,
  lyap: r`\delta(t)\approx\delta_0\,e^{\lambda t}\qquad\Longrightarrow\qquad t_{\text{predict}}\approx\frac{1}{\lambda}\ln\frac{\Delta}{\delta_0}`,
  logh: r`\frac{dt}{ds}=\frac{1}{T+B},\qquad \frac{d\mathbf r_i}{ds}=\frac{\mathbf v_i}{T+B},\qquad \frac{d\mathbf v_i}{ds}=\frac{\mathbf a_i}{U},\qquad
B=-E_0,\quad U=\sum_{i<j}\frac{G m_i m_j}{r_{ij}}`,
};
const IEQ = {
  ri: r`\mathbf r_i`,
  mi: r`m_i`,
  sides: r`\mathbf s_1=\mathbf r_3-\mathbf r_2,\ \mathbf s_2=\mathbf r_1-\mathbf r_3,\ \mathbf s_3=\mathbf r_2-\mathbf r_1`,
};

function tex(src, display) {
  const node = doc.convert(src, { display });
  let html = adaptor.outerHTML(node);
  const label = src.replace(/\s+/g, " ").replace(/"/g, "&quot;");
  return html.replace("<mjx-container", `<mjx-container role="math" aria-label="${label}"`);
}

const here = __dirname;
let page = fs.readFileSync(path.join(here, "src/page.html"), "utf8");
const app = fs.readFileSync(path.join(here, "src/app.js"), "utf8");
const solver = fs.readFileSync(path.join(here, "solver.js"), "utf8");
const orbits = fs.readFileSync(path.join(here, "../data/periodic_orbits.json"), "utf8");
const fatePath = path.join(here, "fatemap.json");
const fate = fs.existsSync(fatePath) ? fs.readFileSync(fatePath, "utf8") : "null";

page = page.replace(/\{\{EQ:(\w+)\}\}/g, (_, k) => { if (!EQ[k]) throw new Error("missing EQ " + k); return tex(EQ[k], true); });
page = page.replace(/\{\{IEQ:(\w+)\}\}/g, (_, k) => { if (!IEQ[k]) throw new Error("missing IEQ " + k); return tex(IEQ[k], false); });
const mjxCss = `mjx-container { display: inline-block; } mjx-container[display="true"] { display: block; margin: 0; } mjx-container svg { overflow: visible; }`;
page = page.split("</style>").join(mjxCss + "\n</style>");
page = page.split("/*__SOLVER__*/").join(solver);
page = page.split("/*__ORBITS__*/null").join(orbits);
page = page.split("/*__FATE__*/null").join(fate);
page = page.split("/*__APP__*/").join(app);

fs.mkdirSync(path.join(here, "dist"), { recursive: true });
const out = path.join(here, "dist/index.html");
fs.writeFileSync(out, page);
console.log(out, (page.length / 1024).toFixed(0) + " KB");
