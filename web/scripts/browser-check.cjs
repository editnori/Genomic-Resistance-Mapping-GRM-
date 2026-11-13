// Lane O1 browser check for the ADS upgrade, run on Vast against a built
// artifact that a vite preview serves with GRM_REPORT_DIR set to one
// read-only fixture.
// node browser-check-final.cjs <playwright-core> <output-dir> <base-url> <fixture-dir>
// GRM_CHROMIUM, when set, is the Chromium executable to launch; unset or
// empty, playwright-core launches its own browser. The check used to pin
// the headless shell of Chromium build 1248 on the release host; to keep that
// build there, set
// GRM_CHROMIUM=/root/.cache/ms-playwright/chromium_headless_shell-1248/chrome-headless-shell-linux64/chrome-headless-shell
// Expectations come from the fixture's own index.json, not from the UI code.
// Assertions that encoded r13's removals are rewritten to the DESIGN.md
// contract: visible titles, the view nav, judged cells as plain text, the
// sort glyph slot, the status line, blue headers, the 1560px rail. Controls
// are checked for names, current state and hit size, not for a pill shape,
// which the user rejected; the primer and drug-class marks have their own
// checks. So do Present's whole deck title and outline label, its footer rows
// and fixed-figure gap, the primer's split tracks, and every species'
// population map (compact below seven blocks).
const fs = require('node:fs');
const path = require('node:path');
const { chromium } = require(process.argv[2]);
const OUT = process.argv[3], BASE = process.argv[4], FIXTURE = process.argv[5];
fs.mkdirSync(OUT, { recursive: true });
const manifest = JSON.parse(fs.readFileSync(path.join(FIXTURE, 'index.json')));
const checks = [], screenshots = [], errors = [];
const check = (name, pass, actual) => checks.push({ name, pass: Boolean(pass), actual });
// Checks a tree gives nothing to test, named in the receipt rather than
// dropped from the total without a word.
const skipped = [];
const median = values => { const s = [...values].sort((a, b) => a - b), k = s.length >> 1; return s.length % 2 ? s[k] : (s[k - 1] + s[k]) / 2; };
const best = d => d.learners.find(l => l.learner === d.best_learner)?.random;
const expected = {
  datasets: manifest.datasets.length,
  paired: manifest.datasets.filter(d => d.best_random !== null && d.best_lineage !== null),
  catalog: manifest.datasets.filter(d => d.catalog !== null && d.best_lineage !== null),
  mic: manifest.datasets.filter(d => d.mic_essential_agreement !== null),
};
const byId = new Map(manifest.datasets.map(d => [d.id, d]));
const resistantOf = p => (p.counts ? p.counts.tp + p.counts.fn : (byId.get(p.dataset)?.n_resistant ?? 0)) + (p.no_call?.resistant ?? 0);
const panelIds = [...new Set(manifest.panels.map(p => p.panel))];
const panelSize = id => manifest.panels.filter(p => p.panel === id).reduce((a, p) => a + resistantOf(p), 0);
const leadId = panelIds.includes('carbapenemase') ? 'carbapenemase' : panelIds.sort((a, b) => panelSize(b) - panelSize(a))[0];
expected.panels = manifest.panels.filter(p => p.panel === leadId);
const gaps = expected.paired.map(d => d.best_random - d.best_lineage);
expected.stats = {
  paired: expected.paired.length,
  medianRandom: median(expected.paired.map(d => d.best_random)),
  medianLineage: median(expected.paired.map(d => d.best_lineage)),
  medianGap: median(gaps),
  above5: gaps.filter(v => v > 0.05).length,
  catalogAhead: expected.catalog.filter(d => d.catalog.balanced_accuracy - d.best_lineage > 0.02).length,
  catalogWithin: expected.catalog.filter(d => Math.abs(d.catalog.balanced_accuracy - d.best_lineage) <= 0.02).length,
  catalogBehind: expected.catalog.filter(d => d.catalog.balanced_accuracy - d.best_lineage < -0.02).length,
  micAbove90: expected.mic.filter(d => d.mic_essential_agreement >= 0.9).length,
  vmeMeets: manifest.datasets.filter(d => best(d)?.vme <= 0.015).length,
  meMeets: manifest.datasets.filter(d => best(d)?.me <= 0.03).length,
  leadPanel: leadId,
  leadPanelDatasets: expected.panels.length,
  // Datasets whose best learner has both error rates with lineages held out
  // and called at least one genome of each class there.
  barsCounted: manifest.datasets.filter(d => { const r = d.learners.find(l => l.learner === d.best_learner)?.lineage; return r && r.vme != null && r.me != null && d.n_resistant - (r.no_call?.resistant ?? 0) > 0 && d.n - d.n_resistant - (r.no_call?.susceptible ?? 0) > 0; }).length,
};
const SHORT = [[/lacks one or both classes/, 'Fold lacks a class'], [/fewer than twice as many lineages/, 'Too few lineages']];
const stateRows = [
  ...SHORT.flatMap(([re, word]) => { const d = manifest.datasets.find(x => x.best_lineage === null && x.best_learner !== null && re.test(x.lineage_unavailable_reason ?? '')); return d ? [[d, word]] : []; }),
  ...(() => { const d = manifest.datasets.find(x => x.catalog === null); return d ? [[d, 'No rule']] : []; })(),
  ...(() => { const d = manifest.datasets.find(x => x.mic_essential_agreement === null); return d ? [[d, 'No MIC']] : []; })(),
];
// A lone rule of 80 characters or more, longer than any rule column at the
// check's widths holds on one line (--shape full writes one).
const longRule = manifest.datasets.find(d => d.scm_rules.length === 1 && d.scm_rules[0].length >= 80);
// The names a species may go by, from the fixture's own words: in full, or
// as its genus initial joined to the rest by a non-breaking space, so the
// short form never ends a line ("Synthetic species C", "S. species C").
const speciesNames = species => { const [genus, ...rest] = species.split(' '); return rest.length > 0 ? [species, genus.charAt(0) + '.\u00a0' + rest.join(' ')] : [species]; };
// A dataset as the app titles it (datasetTitle in format.ts): the species'
// everyday name where it has one, else its short name, then "vs" held to
// the drug by a non-breaking space ("S. species C vs ceftriaxone").
const EVERYDAY = new Map([['Mycobacterium tuberculosis', 'TB'], ['Neisseria gonorrhoeae', 'Gonorrhea'], ['Salmonella enterica', 'Salmonella'], ['Staphylococcus aureus', 'Staph'], ['Streptococcus pneumoniae', 'Pneumococcus']]);
const datasetTitle = d => (EVERYDAY.get(d.species) ?? speciesNames(d.species).at(-1)) + ' vs\u00a0' + d.antibiotic;
// The dataset and error view whose interval text runs widest, in characters
// as the learner table prints it (one decimal), from the tree's own
// intervals/; null for a tree without them.
const widestInterval = (() => {
  const dir = path.join(FIXTURE, 'intervals');
  if (!fs.existsSync(dir)) return null;
  const text = ci => (ci ? (ci.lo * 100).toFixed(1) + '%–' + (ci.hi * 100).toFixed(1) + '%' : '');
  let widest = null;
  for (const name of fs.readdirSync(dir).sort()) {
    const id = name.replace(/\.json$/, '');
    if (!byId.has(id)) continue;
    for (const s of JSON.parse(fs.readFileSync(path.join(dir, name))).splits) {
      const chars = Math.max(0, ...s.learners.map(l => text(l.vme).length)) + Math.max(0, ...s.learners.map(l => text(l.me).length));
      if (widest === null || chars > widest.chars) widest = { id, scheme: s.scheme, chars };
    }
  }
  return widest;
})();
// From the tree's own datasets/: the dataset whose decision tree, at the
// inspector's default three levels, has the most leaf columns, and the one
// whose catalog gene groups hold the longest class or subclass list (one
// slash-joined run, such as "lincosamide/oxazolidinone/..."); null when no
// dataset has one. Also the trees whose labels the label check reads: the
// widest, then those with the most default-view split labels longer than 18
// characters (the app's own words for a feature, featureLabel's: catalog
// families, else gene, else product), then the most splits, by id; LABEL_TREES in all.
const LABEL_TREES = 6;
const widestPanels = (() => {
  let tree = null, catalog = null;
  const ranked = [];
  for (const d of manifest.datasets) {
    const file = path.join(FIXTURE, 'datasets', d.id + '.json');
    if (!fs.existsSync(file)) continue;
    const r = JSON.parse(fs.readFileSync(file));
    const nodes = r.models?.cart?.nodes ?? [], depth = Math.min(3, r.models?.cart_depth ?? 0);
    const cards = new Map((r.features ?? []).map(f => [f.feature, f]));
    const families = c => { const names = [...new Set(c?.catalog_families ?? [])].sort(); return names.length === 0 ? null : names.length === 1 ? names[0] : `k-mers matching ${names.length} catalog families (${names.join(', ')})`; };
    const words = i => { const c = cards.get(i); return families(c) ?? c?.locus?.gene ?? c?.locus?.product ?? 'k-mer feature ' + i; };
    let columns = 0, splits = 0, long = 0;
    const visit = i => { const n = nodes[i]; if (n === undefined) return; if (n.split === null || n.depth >= depth) { columns += 1; return; } splits += 1; if (words(n.split.feature).length > 18) long += 1; visit(n.split.absent); visit(n.split.present); };
    if (nodes.length > 0) visit(0);
    if (columns > 0 && (tree === null || columns > tree.columns)) tree = { id: d.id, columns };
    if (splits > 0) ranked.push({ id: d.id, long, splits });
    const chars = Math.max(0, ...(r.catalog?.groups ?? []).map(g => Math.max(g.class.length, g.subclass.length)));
    if (chars > 0 && (catalog === null || chars > catalog.chars)) catalog = { id: d.id, chars, groups: r.catalog.groups.length };
  }
  ranked.sort((a, b) => b.long - a.long || b.splits - a.splits || (a.id < b.id ? -1 : 1));
  const labels = tree === null ? [] : [tree.id, ...ranked.map(x => x.id).filter(id => id !== tree.id)].slice(0, LABEL_TREES);
  return { tree, catalog, labels };
})();
const inspectId = manifest.datasets.some(d => d.id === 'escherichia_coli__ceftriaxone') ? 'escherichia_coli__ceftriaxone' : manifest.datasets[0].id;
const speciesSlug = manifest.species.find(s => s.slug === 'escherichia_coli' && s.genomes_analyzed > 0)?.slug ?? manifest.species.find(s => s.genomes_analyzed > 0).slug;
const speciesName = manifest.species.find(s => s.slug === speciesSlug).species;
const allValid = manifest.datasets.every(d => d.evaluation?.valid === true && (manifest.evaluation_revision === undefined || d.evaluation.revision === manifest.evaluation_revision));
const CHIP = allValid ? 'Tested within one public database only' : 'Provisional metrics';
const isLimitsHref = href => href === '/#limits' || href === '#limits' || href === '#m-limits';
// How long a check waits for the state it reads (a transition that lands, a
// strip that redraws, a link that follows the scroll) before it reads what is
// there and judges that. A wait on state returns as soon as the state holds,
// so a fast machine reads as early as before; the deadline only decides how
// slow a machine may be before the check calls the state missing.
const STATE_MS = 10000;
// The inspector's freshness window for a launch point (FRESH_MS in
// web/src/ui/launch.ts), read from the app's own source so the two cannot
// drift apart.
const FRESH_MS = Number(fs.readFileSync(path.join(__dirname, '../src/ui/launch.ts'), 'utf8').match(/const FRESH_MS = (\d+);/)[1]);
// When the title slide's last beat lands: the sum of its stage delays
// (useSlideTimeline's delays for the title in web/src/present/PresentMode.tsx).
const TITLE_REVEAL_MS = fs.readFileSync(path.join(__dirname, '../src/present/PresentMode.tsx'), 'utf8').match(/slide\.kind === "title" \? \[([\d, ]+)\]/)[1].split(',').reduce((a, v) => a + Number(v), 0);
const SLIDE_TITLES = ['Genomic Resistance Mapping', 'Corpus', 'Lineage cost', 'Known-gene catalog', null, 'MIC prediction', 'One interpretable rule'];

const settle = async page => {
  await page.evaluate(() => document.fonts.ready);
  await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
};
const shot = async (page, name, full = false) => {
  await settle(page);
  const file = path.join(OUT, name + '.png');
  await page.screenshot({ path: file, fullPage: full });
  screenshots.push(file);
  console.log('Screenshot ' + name);
};
const tallShot = async (page, name, width, height) => {
  await settle(page);
  const h = await page.evaluate(() => { const v = document.querySelector('[data-view-scroll]'); return v ? v.scrollHeight + v.getBoundingClientRect().top + 16 : document.documentElement.scrollHeight; });
  const cdp = await page.context().newCDPSession(page);
  await cdp.send('Emulation.setDeviceMetricsOverride', { width, height: Math.min(Math.ceil(h), 16000), deviceScaleFactor: 1, mobile: false });
  await page.setViewportSize({ width, height: Math.min(Math.ceil(h), 16000) });
  await shot(page, name);
  await page.setViewportSize({ width, height });
  await cdp.send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 1, mobile: false });
  await cdp.detach();
};

const measure = () => {
  const main = document.querySelector('dialog.presentation-dialog main');
  const section = main?.querySelector('section');
  const rect = e => e ? ({ x: e.getBoundingClientRect().x, y: e.getBoundingClientRect().y, width: e.clientWidth, height: e.clientHeight, scrollWidth: e.scrollWidth, scrollHeight: e.scrollHeight, bottom: e.getBoundingClientRect().bottom }) : null;
  const texts = [...(section?.querySelectorAll('svg text') ?? [])].map(e => ({ text: e.textContent, box: e.getBoundingClientRect().toJSON(), svg: e.closest('svg').getBoundingClientRect().toJSON() }));
  const textOutside = texts.filter(t => t.box.left < t.svg.left - 1 || t.box.right > t.svg.right + 1 || t.box.top < t.svg.top - 1 || t.box.bottom > t.svg.bottom + 1);
  const caption = section?.querySelector('[data-caption]');
  const figure = section?.querySelector('[data-figure]');
  const fit = section?.querySelector('[data-figure-fit]');
  const visible = e => !e.closest('[role="tooltip"]') && !e.classList.contains('sr-only') && e.getClientRects().length > 0;
  const contentBottom = Math.max(0, ...[...(figure?.querySelectorAll('*') ?? [])].filter(visible).map(e => e.getBoundingClientRect().bottom));
  const figureBox = figure?.getBoundingClientRect();
  return { viewport: [innerWidth, innerHeight], main: rect(main), section: rect(section), fit: rect(fit), contentBottom, caption: rect(caption), figureTop: figureBox?.top ?? null, mainTop: main?.getBoundingClientRect().top ?? null, textOutside, marks: section?.querySelectorAll('[data-dataset]').length ?? 0, titleX: section?.querySelector('h2')?.getBoundingClientRect().x ?? null, brandX: document.querySelector('dialog.presentation-dialog header svg')?.getBoundingClientRect().x ?? null };
};

async function referenceChecks(page, tag, finding, panelSeen = [], panelPaged = false) {
  if (![1, 2, 3, 4].includes(finding)) return;
  const geometry = await page.locator('dialog.presentation-dialog main svg:not(.grm-mark)').first().evaluate(svg => ({ width: svg.getBoundingClientRect().width, height: svg.getBoundingClientRect().height, refs: [...svg.querySelectorAll('line[aria-hidden="true"]')].map(l => ({ x: Number(l.getAttribute('x1')), y: Number(l.getAttribute('y1')) })), ids: [...svg.querySelectorAll('[data-dataset]')].map(e => e.getAttribute('data-dataset')), text: svg.textContent }));
  const close = (actual, want) => Math.abs(actual - want) <= 1.5;
  let wanted = [];
  if (finding === 1) {
    const lo = Math.floor(Math.min(-0.02, ...gaps) * 100) / 100, hi = Math.ceil(Math.max(0.1, ...gaps) * 100) / 100;
    wanted = [0.05, 0].map(v => 20 + (v - lo) / (hi - lo) * (geometry.width - 40));
    check(tag + ' lineage region counts', new RegExp((gaps.length - expected.stats.above5) + ' ≤5( points)?').test(geometry.text) && new RegExp(expected.stats.above5 + ' >5( points)?').test(geometry.text), geometry.text);
  }
  if (finding === 4) {
    const values = expected.mic.map(d => d.mic_essential_agreement), lo = Math.min(0.8, Math.floor(Math.min(...values) * 10) / 10);
    wanted = [0.9].map(v => 20 + (v - lo) / (1 - lo) * (geometry.width - 40));
    check(tag + ' MIC region counts', geometry.text.includes((values.length - expected.stats.micAbove90) + ' <90%') && geometry.text.includes(expected.stats.micAbove90 + ' ≥90%'), geometry.text);
  }
  if (finding === 3) {
    // A report from before revision 11 has no panel counts. The page then
    // rebuilds them from each rate and the dataset's class sizes (joinPanels).
    const counts = expected.panels.map(p => {
      const called = p.counts ? p.counts.tp + p.counts.fn : (byId.get(p.dataset)?.n_resistant ?? 0), resistant = called + (p.no_call?.resistant ?? 0), m = p.best_model_counts;
      return { caught: p.counts ? p.counts.tp : Math.round(p.metrics.sensitivity * called), called, modelCaught: m ? m.tp : Math.round(p.best_model.sensitivity * resistant), modelResistant: m ? m.tp + m.fn : resistant };
    });
    const panel = counts.reduce((a, c) => a + c.caught, 0) / counts.reduce((a, c) => a + c.called, 0);
    const model = counts.reduce((a, c) => a + c.modelCaught, 0) / counts.reduce((a, c) => a + c.modelResistant, 0);
    const start = Math.min(220, Math.round(geometry.width * 0.5)) + 10;
    wanted = [panel, model].map(v => start + v * (geometry.width - 20 - start));
    // Pooled reference lines span every row on every page; a page draws
    // only lead-panel datasets, and the pages together are checked in present().
    const lead = new Set(expected.panels.map(p => p.dataset));
    if (panelPaged) check(tag + ' lead panel datasets on this page', geometry.ids.length > 0 && geometry.ids.every(id => lead.has(id) && !panelSeen.includes(id)), geometry.ids);
    else check(tag + ' lead panel datasets', JSON.stringify([...geometry.ids].sort()) === JSON.stringify(expected.panels.map(p => p.dataset).sort()), geometry.ids);
    panelSeen.push(...geometry.ids);
  }
  if (finding === 2) {
    const lo = Math.max(0, Math.floor(Math.min(0.5, ...expected.catalog.flatMap(d => [d.catalog.balanced_accuracy, d.best_lineage])) * 10) / 10);
    // CatalogScatter draws the 50% line only when its axis starts below 50%
    // (lo < 0.5); at lo = 0.5 the axis floor is 50% and no line is drawn.
    // A missing line counts only from a figure that drew: a box, and one
    // mark per dataset with a catalog score and a lineage score.
    const at50 = geometry.refs.filter(r => close(r.y, 22 + (1 - 0.5) / (1 - lo) * (geometry.height - 56)));
    const drawn = geometry.width > 0 && geometry.height > 0 && geometry.ids.length === expected.catalog.length;
    // RefLineH is the figure's only aria-hidden line, so with no 50% line
    // there is no reference line at all.
    check(tag + ' catalog 50% reference', drawn && (lo < 0.5 ? at50.length > 0 : geometry.refs.length === 0), { lo, drawn, marks: geometry.ids.length, refs: geometry.refs });
    check(tag + ' catalog region counts', geometry.text.includes(expected.stats.catalogAhead + ' ahead by more than 2') && geometry.text.includes(expected.stats.catalogWithin + ' within 2 points') && geometry.text.includes(expected.stats.catalogBehind + ' behind by more than 2'), geometry.text);
  } else check(tag + ' finding ' + finding + ' reference positions', wanted.length === geometry.refs.length && wanted.every((x, i) => close(geometry.refs[i].x, x)), { wanted, actual: geometry.refs });
}

/** Judged values are plain text with accessible state words, not pills. */
async function judgedCheck(page, tag) {
  const values = await page.locator('[data-verdict]').evaluateAll(cells => {
    const rgb = s => (s.match(/[\d.]+/g) ?? []).map(Number);
    const lum = c => c.slice(0, 3).map(x => { const v = x / 255; return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4; }).reduce((s, v, i) => s + v * [0.2126, 0.7152, 0.0722][i], 0);
    const painted = e => {
      const layers = [];
      for (let n = e; n; n = n.parentElement) { const c = rgb(getComputedStyle(n).backgroundColor); if (c.length >= 3) layers.push(c); }
      const paper = rgb(getComputedStyle(document.documentElement).getPropertyValue('--paper'));
      let under = paper.length >= 3 ? paper.slice(0, 3) : [255, 255, 255];
      for (const c of layers.reverse()) { const a = c[3] ?? 1; under = under.map((v, i) => c[i] * a + v * (1 - a)); }
      return under;
    };
    return cells.map(e => {
      const s = getComputedStyle(e), a = lum(rgb(s.color)), b = lum(painted(e));
      return { text: e.textContent.trim(), verdict: e.getAttribute('data-verdict'), border: s.borderTopWidth, shadow: s.boxShadow, radius: parseFloat(s.borderTopLeftRadius), ratio: (Math.max(a, b) + 0.05) / (Math.min(a, b) + 0.05) };
    });
  });
  check(tag + ' judged values are plain text with 4.5:1 contrast', values.length > 0 && values.every(v => v.border === '0px' && v.shadow === 'none' && v.radius === 0 && v.ratio >= 4.5), values);
}

/** Top-bar and page controls: each has a name and at least a 24px target,
 *  and none is drawn as a pill (radius at half its height or more). */
async function controlCheck(page, tag) {
  const controls = await page.evaluate(() => [...document.querySelectorAll('header button, header [role="combobox"], [data-view-scroll] button.ui-control-target, [data-view-scroll] [role="combobox"]')]
    .filter(e => e.getClientRects().length > 0 && !e.closest('table') && !e.matches('[data-row-open]') && !e.closest('[data-view-options-panel] > div:first-child'))
    .map(e => ({ name: (e.getAttribute('aria-label') ?? e.textContent ?? '').trim().slice(0, 30), r: parseFloat(getComputedStyle(e).borderTopLeftRadius), h: e.getBoundingClientRect().height, w: e.getBoundingClientRect().width })));
  check(tag + ' controls are named with 24px targets', controls.length > 0 && controls.every(c => c.name.length > 0 && c.h >= 24 - 0.5 && c.w >= 24 - 0.5), controls.filter(c => !(c.name.length > 0 && c.h >= 23.5 && c.w >= 23.5)));
  check(tag + ' no control is a pill', controls.every(c => c.r < c.h / 2 - 1), controls.filter(c => c.r >= c.h / 2 - 1));
}

/** The brand mark and the view's title start on one left edge (the 24px
 *  gutter both share) where the 1560px rail is not centred. */
async function edgeCheck(page, tag, selector) {
  const edges = await page.evaluate(sel => ({ brand: document.querySelector('header a[href="/"] svg')?.getBoundingClientRect().x, first: document.querySelector(sel)?.getBoundingClientRect().x, width: innerWidth }), selector);
  const applies = edges.width <= 1560 + 16 + 48;
  check(tag + ' left edge ' + selector, !applies || (edges.brand !== undefined && edges.first !== undefined && Math.abs(edges.brand - edges.first) <= 1), edges);
}

/** The Overview primer: five steps above the numbered findings,
 *  a real report as its example, and a skip link that lands focus on the
 *  first finding's heading. */
async function primerCheck(page, tag) {
  const p = await page.evaluate(() => {
    const primer = document.getElementById('primer');
    const corpus = document.getElementById('corpus');
    const example = primer ? [...primer.querySelectorAll('a')].find(a => !a.getAttribute('href')?.startsWith('#')) : null;
    return {
      steps: primer ? [...primer.querySelectorAll('h3')].map(h => h.textContent) : [],
      before: primer && corpus ? Boolean(primer.compareDocumentPosition(corpus) & Node.DOCUMENT_POSITION_FOLLOWING) : false,
      figures: primer ? [...primer.querySelectorAll('[data-primer-figure]')].map(f => ({ id: f.getAttribute('data-primer-figure'), right: f.getBoundingClientRect().right, view: document.querySelector('[data-view-scroll]').getBoundingClientRect().right })) : [],
      example: example?.getAttribute('href') ?? null,
    };
  });
  check(tag + ' primer: five steps above the findings', p.before && JSON.stringify(p.steps) === JSON.stringify(['A dataset', 'A genome feature', 'Two tests', 'Reading one result', 'What the numbers support']), p);
  check(tag + ' primer figures stay inside the view', p.figures.length >= 3 && p.figures.every(f => f.right <= f.view + 1), p.figures);
  check(tag + ' primer example opens a real report', p.example !== null && manifest.datasets.some(d => p.example.includes('dataset=' + d.id)), p.example);
  await splitTracks(page, tag + ' primer', '#primer');
  await page.getByRole('link', { name: /Skip to the findings/ }).click();
  await page.waitForFunction(() => document.activeElement?.id === 'corpus-title');
  check(tag + ' primer skip link moves focus to the first finding', await page.evaluate(() => document.activeElement?.id) === 'corpus-title', await activeName(page));
}

/** The "Two tests" drawing: in each split its name and score share a line,
 *  its five lineage groups sit on one line, and the two splits' groups start
 *  at the same x, so the rows read group for group. */
async function splitTracks(page, tag, scope) {
  const t = await page.evaluate(sel => {
    const root = document.querySelector(sel);
    return [...(root?.querySelectorAll('[data-primer-figure="splits"] [data-split-row]') ?? [])].map(row => {
      const head = row.firstElementChild;
      const [label, score] = [...head.children].map(e => e.getBoundingClientRect());
      return { groups: [...row.querySelectorAll('[data-lineage-group]')].map(g => ({ left: Math.round(g.getBoundingClientRect().left * 10) / 10, top: Math.round(g.getBoundingClientRect().top * 10) / 10 })), labelTop: label?.top, scoreTop: score?.top, scoreBottom: score?.bottom, labelBottom: label?.bottom, figureRight: row.closest('figure').getBoundingClientRect().right, scoreRight: score?.right };
    });
  }, scope);
  const oneLine = t.length === 2 && t.every(r => r.groups.length === 5 && r.groups.every(g => Math.abs(g.top - r.groups[0].top) <= 0.5));
  const aligned = oneLine && t[0].groups.every((g, i) => Math.abs(g.left - t[1].groups[i].left) <= 0.5);
  const heads = t.every(r => r.scoreTop !== undefined && r.scoreTop < r.labelBottom && r.labelTop < r.scoreBottom && r.scoreRight <= r.figureRight + 0.5);
  check(tag + ' two-tests squares on even tracks', oneLine && aligned && heads, t);
}

/** Every analyzed species' population map: from COMPACT_BLOCKS (7) blocks
 *  up the columns spread; below it the map is compact, columns of 24 to
 *  72px, with a visible note and one table row per block. */
async function smallMaps(page, tag) {
  const compact = [];
  for (const s of manifest.species.filter(x => x.genomes_analyzed > 0)) {
    await page.goto(BASE + '/species/' + s.slug);
    await page.getByRole('grid', { name: /lineage block/ }).first().waitFor();
    await settle(page);
    const m = await page.evaluate(() => {
      const map = document.querySelector('#sp-structure [data-block-layout]');
      const grid = map?.querySelector('[role="grid"]');
      const blocks = Number(((grid?.getAttribute('aria-label') ?? '').replace(/,/g, '').match(/^(\d+) lineage blocks? of/) ?? [])[1] ?? NaN);
      const box = map?.querySelector('[data-compact-blocks]');
      const plot = grid?.parentElement, scroller = plot?.parentElement;
      const note = box?.querySelector('p');
      const wrap = box?.querySelector('table')?.parentElement;
      return {
        layout: map?.getAttribute('data-block-layout') ?? null, blocks,
        compact: box ? Number(box.getAttribute('data-compact-blocks')) : null,
        rows: box ? box.querySelectorAll('tbody tr').length : 0,
        note: note?.textContent ?? null, noteShown: note ? note.getBoundingClientRect().height > 0 : false,
        plotW: plot?.getBoundingClientRect().width ?? null, viewW: scroller?.clientWidth ?? null,
        tableFits: wrap ? wrap.scrollWidth <= wrap.clientWidth + 1 : null,
        scrollNote: [...map.querySelectorAll('p')].some(p => /scroll sideways for the rest/.test(p.textContent) && p.getBoundingClientRect().height > 0),
        doc: document.documentElement.scrollWidth - innerWidth,
      };
    });
    // Compact columns are 24 to 72px: they fit the view when it holds 24px
    // per block, and otherwise scroll with the visible sideways note, as a
    // spread map does.
    const small = m.blocks >= 1 && m.blocks < 7;
    const width = m.plotW >= 24 * m.blocks - 0.5 && m.plotW <= 72 * m.blocks + 0.5 && (m.viewW >= 24 * m.blocks ? m.plotW <= m.viewW + 1 : m.scrollNote);
    const ok = m.blocks >= 1 && m.doc <= 1 && (small
      ? m.layout === 'compact' && m.compact === m.blocks && m.rows === m.blocks && m.noteShown && m.note.includes(`make only ${m.blocks} block`) && width && m.tableFits === true
      : m.layout === 'spread' && m.compact === null);
    check(tag + ' species ' + s.slug + ' ' + (small ? 'compact' : 'spread') + ' population map', ok, { slug: s.slug, ...m });
    if (small) {
      compact.push(s.slug + ':' + m.blocks);
      await page.locator('#sp-structure').evaluate(e => e.scrollIntoView({ block: 'start' }));
      await shot(page, tag + '-small-map-' + s.slug);
    }
  }
  check(tag + ' small maps listed', true, compact);
}

async function reachCheck(page, tag) {
  const over = await page.evaluate(() => document.documentElement.scrollWidth > innerWidth + 1 ? [{ width: innerWidth, scroll: document.documentElement.scrollWidth }] : []);
  check(tag + ' nothing reaches past the page edge', over.length === 0, over);
}

/** Every view opens with a visible 32px title and one quiet status line
 *  that carries the audit chip and the report's date. */
async function headerCheck(page, tag, title) {
  const h = await page.evaluate(() => {
    const h1 = document.querySelector('[data-view-scroll] h1');
    const header = h1?.closest('header');
    return { text: h1?.textContent, size: h1 ? parseFloat(getComputedStyle(h1).fontSize) : null, srOnly: h1?.classList.contains('sr-only') ?? null, status: header?.textContent ?? '', scopeHref: header?.querySelector('a[data-audit-status], a[href$="#limits"], a[href="#m-limits"]')?.getAttribute('href') ?? null };
  });
  check(tag + ' visible title', h.text === title && h.size === 32 && h.srOnly === false, h);
  // On a provisional report, Overview and Methods print the status in bold with the reasons
  // under the line, and carry no scope link (ReportLine in findings.tsx).
  const bare = !allValid && / (overview|methods)$/.test(tag);
  check(tag + ' plain status line with scope and report date', h.status.includes(CHIP) && h.status.includes('Report built') && !h.status.includes('BV-BRC lab records') && (bare ? h.scopeHref === null : isLimitsHref(h.scopeHref)), h);
}

// The guidance's own acceptance wording, read from the release summary
// whose revision, generation time and dataset count match this tree:
// reference/release-summary/<release>/guidance.json, written by
// scripts/release-summary/guidance.py. A synthetic tree matches none, and
// the checks that compare with it are skipped by name.
const GUIDANCE = (() => {
  const dir = path.join(__dirname, '../../reference/release-summary');
  const names = fs.existsSync(dir) ? fs.readdirSync(dir).sort() : [];
  for (const name of names) {
    const file = path.join(dir, name, 'guidance.json');
    if (!fs.existsSync(file)) continue;
    const g = JSON.parse(fs.readFileSync(file));
    if (g.revision === manifest.evaluation_revision && g.generated_at === manifest.generated_at && g.datasets === manifest.datasets.length) return { release: name, rows: g.rows };
  }
  return null;
})();
const GUIDANCE_MISSING = 'no reference/release-summary/*/guidance.json has this tree\'s revision ' + manifest.evaluation_revision + ', generation time ' + manifest.generated_at + ' and ' + manifest.datasets.length + ' datasets';
/** One row of the release summary as the page counts it: datasets with both
 *  rates, those that meet the wording, and the M. tuberculosis datasets
 *  among those. */
const guidanceRow = key => {
  const r = GUIDANCE?.rows?.[key];
  return r === undefined ? null : { n: r.n, meets: r.meets_wording, tuberculosis: r.meets_wording - r.without_tuberculosis.meets_wording };
};
const WORDING_WHO = { 'Best learner for each dataset': 'best', 'Set covering machine': 'scm', 'Decision tree (CART)': 'cart', 'Random forest': 'rf', 'Gradient boosted trees': 'gbt', 'L1 logistic regression': 'lasso' };
const WORDING_SPLIT = { 'Lineages held out': 'lineage', 'Random split': 'random' };
const toCount = s => (s === undefined || s === null ? null : Number(String(s).replace(/,/g, '')));

/** Overview's bars section, its table open: the last column of "Each
 *  learner against the bars" and its Datasets column, for each of the
 *  twelve rows, and the counts in the note on the guidance's wording, all
 *  equal to the release summary's rows["<who>:<split>"].meets_wording and
 *  .n. */
async function wordingCheck(page, tag) {
  const tableName = tag + ' overview counts by the guidance\'s wording match the release summary, each learner and split';
  const noteName = tag + ' overview note on the guidance\'s wording gives the release summary\'s counts';
  if (GUIDANCE === null) { for (const name of [tableName, noteName]) skipped.push({ check: name, reason: GUIDANCE_MISSING }); return; }
  const shown = await page.evaluate(() => {
    const region = document.querySelector('#bars [role="region"][aria-label="Each learner against the bars"]');
    const head = region === null ? [] : [...region.querySelectorAll('thead th')].map(t => t.textContent.trim());
    const col = head.indexOf('Meets the guidance\'s wording');
    const rows = [];
    let group = null;
    for (const tr of region === null ? [] : region.querySelectorAll('tbody tr')) {
      const g = tr.querySelector('th[scope="rowgroup"]');
      if (g !== null) { group = g.textContent.trim(); continue; }
      const cells = [...tr.children].map(c => c.textContent.trim());
      rows.push({ group, who: cells[0], n: cells[1], meets: col < 0 ? null : cells[col] });
    }
    const note = [...document.querySelectorAll('#bars p')].map(p => p.textContent.replace(/\s+/g, ' ').trim()).find(t => t.startsWith('The guidance\'s own acceptance wording')) ?? null;
    return { head, col, rows, note };
  });
  const rows = shown.rows.map(r => {
    const key = WORDING_WHO[r.who] + ':' + WORDING_SPLIT[r.group];
    const want = guidanceRow(key);
    return { key, n: toCount(r.n), meets: toCount(r.meets), want: want === null ? null : { n: want.n, meets: want.meets } };
  });
  check(tableName, shown.col > 0 && rows.length === 12 && new Set(rows.map(r => r.key)).size === 12 && rows.every(r => r.want !== null && r.n === r.want.n && r.meets === r.want.meets), { release: GUIDANCE.release, column: shown.col, rows });
  const lineage = guidanceRow('best:lineage'), random = guidanceRow('best:random');
  const t = shown.note ?? '';
  const l = t.match(/meets on ([\d,]+) of ([\d,]+) datasets? with lineages held out/), r = t.match(/(?:meets|and) on ([\d,]+) of ([\d,]+)(?: datasets?)? on random splits/);
  const tb = /no M\. tuberculosis dataset is among them/.test(t) ? [0, 0] : (m => (m === null ? null : [toCount(m[1]), toCount(m[2])]))(t.match(/Of those, ([\d,]+) with lineages held out and ([\d,]+) on random splits/));
  const read = { lineage: l === null ? null : [toCount(l[1]), toCount(l[2])], random: r === null ? null : [toCount(r[1]), toCount(r[2])], tuberculosis: tb ?? (lineage.meets + random.meets === 0 ? [0, 0] : null) };
  const want = { lineage: [lineage.meets, lineage.n], random: [random.meets, random.n], tuberculosis: [lineage.tuberculosis, random.tuberculosis] };
  check(noteName, shown.note !== null && JSON.stringify(read) === JSON.stringify(want), { release: GUIDANCE.release, read, want, note: shown.note });
}

// The outside-study rows the app bundles (web/src/data/outside-study.json,
// written by scripts/outside-study-results/site_data.py), read here by the
// documents' rule: a row holds when neither reading is "below", and a row
// with a reading outside the three is not read.
const OUTSIDE = JSON.parse(fs.readFileSync(path.join(__dirname, '../src/data/outside-study.json')));
const OUTSIDE_READ = ['below', 'consistent', 'above'];
const outsideExpected = OUTSIDE.populations.map(p => {
  const known = r => OUTSIDE_READ.includes(r.sensitivity) && OUTSIDE_READ.includes(r.specificity);
  const scm = p.rows.filter(r => r.learner === 'scm' && known(r));
  return { id: p.id, read: scm.length, holds: scm.filter(r => r.sensitivity !== 'below' && r.specificity !== 'below').length, genomes: p.species.reduce((a, x) => a + x.frozen, 0), below: { sensitivity: scm.filter(r => r.sensitivity === 'below').length, specificity: scm.filter(r => r.specificity === 'below').length }, above: { sensitivity: scm.filter(r => r.sensitivity === 'above').length, specificity: scm.filter(r => r.specificity === 'above').length }, ids: scm.map(r => r.dataset) };
});

/** Newly deposited genomes: the lead's two counts, one mark per set
 *  covering machine row read in each of the four plots, hollow where it
 *  reads below, each plot's above and below counts on the line above it and
 *  clear of the drawn area, and the populations in panels of their own. A pinned mark offers Open dataset only when the
 *  loaded index has that dataset. */
async function outsideCheck(page, tag) {
  const section = page.locator('#outside');
  await section.waitFor();
  const read = await section.evaluate(s => ({
    lead: s.querySelector('h2 + p')?.textContent ?? '',
    panels: [...s.querySelectorAll('h3')].map(h => h.textContent),
    plots: [...s.querySelectorAll('[role="listbox"]')].map(l => {
      // The plot's counts line sits above its svg, outside the drawn area.
      const plot = l.closest('[data-outside-plot]');
      const counts = plot?.querySelector('[data-outside-counts]');
      const box = el => { const r = el.getBoundingClientRect(); return { left: r.left, right: r.right, top: r.top, bottom: r.bottom }; };
      const c = counts ? box(counts) : null, v = box(l.querySelector('svg'));
      return {
        label: l.getAttribute('aria-label'),
        marks: [...l.querySelectorAll('[data-outside]')].map(g => [g.getAttribute('data-outside'), g.getAttribute('data-reading'), g.querySelector('circle').getAttribute('fill')]),
        counts: counts?.textContent ?? null,
        apart: c !== null && (c.bottom <= v.top || c.top >= v.bottom || c.right <= v.left || c.left >= v.right),
        svgText: [...l.querySelectorAll('svg text')].map(t => t.textContent).filter(t => /above|below/.test(t)),
      };
    }),
  }));
  const [primary, same] = outsideExpected;
  const fmt = n => n.toLocaleString('en-US');
  check(tag + ' outside lead gives both populations apart', read.lead.includes(fmt(primary.genomes) + ' genomes from studies the cohorts do not hold') && read.lead.includes(primary.holds + ' of ' + primary.read + ' datasets') && read.lead.includes(fmt(same.genomes) + ' new genomes') && read.lead.includes(same.holds + ' of ' + same.read), read.lead);
  check(tag + ' outside panels, one per population', read.panels.length === 2 && read.panels[0].endsWith(primary.read + ' datasets') && read.panels[1].endsWith(same.read + ' datasets'), read.panels);
  const want = [[primary, 'sensitivity'], [primary, 'specificity'], [same, 'sensitivity'], [same, 'specificity']];
  check(tag + ' outside plots: one mark per row read, hollow where it reads below, counts as the rows give', read.plots.length === 4 && want.every(([p, m], i) => {
    const plot = read.plots[i];
    return plot.marks.length === p.read && JSON.stringify(plot.marks.map(x => x[0]).sort()) === JSON.stringify([...p.ids].sort())
      && plot.marks.filter(x => x[1] === 'below').length === p.below[m] && plot.marks.every(x => (x[1] === 'below') === (x[2] === 'rgb(var(--paper-pure))'))
      && plot.counts === p.above[m] + ' read above, ' + p.below[m] + ' read below' && plot.svgText.length === 0 && plot.label.startsWith(m.charAt(0).toUpperCase() + m.slice(1));
  }), read.plots.map(p => ({ label: p.label, marks: p.marks.length, counts: p.counts, svgText: p.svgText })));
  check(tag + ' outside counts sit above each plot, clear of its svg', read.plots.length === 4 && read.plots.every(p => p.apart), read.plots.map(p => ({ label: p.label, apart: p.apart })));
  // Pin the first mark by keyboard: Enter on the plot's listbox pins its
  // first mark, and the strip offers Open dataset only for a dataset the
  // index holds.
  const first = page.locator('#outside [role="listbox"]').first();
  await first.focus(); await page.keyboard.press('Home'); await page.keyboard.press('Enter');
  const firstId = await first.evaluate(l => l.querySelector('[data-outside]')?.getAttribute('data-outside') ?? null);
  const strip = await section.evaluate(s => ({ open: [...s.querySelectorAll('button')].some(b => b.textContent === 'Open dataset'), clear: [...s.querySelectorAll('button')].some(b => b.getAttribute('aria-label') === 'Clear the pinned mark') }));
  check(tag + ' outside pinned mark opens its dataset only where the index has it', strip.clear && strip.open === byId.has(firstId), { firstId, inIndex: byId.has(firstId), ...strip });
  await page.keyboard.press('Escape');
}

// The input comparison the app bundles (web/src/data/input-comparison.json,
// written by scripts/input-comparison/site_data.py from the reference
// reading): one row per changed input, then the default's two repeat runs.
const INPUTS = JSON.parse(fs.readFileSync(path.join(__dirname, '../src/data/input-comparison.json')));
const fixed4 = v => (v < 0 ? '−' + (-v).toFixed(4) : v.toFixed(4));
// The page's table keeps every column and its bars from BARS_FROM up; under
// it, the slide's narrow set (Change, Worse, Better, Sign p) with no bars.
const inputsExpected = INPUTS.arms.map(a => ({ id: a.id, twin: a.twin, cells: [String(a.all.worse), String(a.all.better), String(a.all.datasets), fixed4(a.all.p), fixed4(a.all.median)], narrow: [String(a.all.worse), String(a.all.better), fixed4(a.all.p)] }));
// Names a reader must never see: the reading's own ids and terms.
const INTERNAL = /H-1|\barm\b|sentinel|twin|family C|\bK\b|\bJ\b|cap-\d|grr-v1|budget-\d|k-c|k-h|k-p|\bk2[15]\b/;

/** Input settings against the default: every row with the reading's counts,
 *  sign p and median, the changed inputs and the default's repeat runs in
 *  groups of their own, each bar drawn to its counts on one scale, no
 *  internal name, and the by-species and by-learner tables each behind a
 *  control of its own, opened by keyboard. */
async function inputsCheck(page, tag) {
  const section = page.locator('#inputs');
  await section.waitFor();
  const read = await section.evaluate(s => ({
    lead: s.querySelector('h2 + p')?.textContent ?? '',
    head: [...(s.querySelector('tbody[data-input-group]')?.closest('table')?.querySelectorAll('thead th') ?? [])].map(th => th.textContent),
    fits: (() => { const t = s.querySelector('tbody[data-input-group]')?.closest('table'); const f = t?.parentElement; return t && f ? t.getBoundingClientRect().width <= f.clientWidth + 0.5 : false; })(),
    groups: [...s.querySelectorAll('tbody[data-input-group]')].map(g => ({ title: g.getAttribute('data-input-group'), rows: [...g.querySelectorAll('tr[data-input]')].map(tr => ({ id: tr.getAttribute('data-input'), cells: [...tr.querySelectorAll('td')].map(td => td.textContent).filter(t => t !== ''), bars: [...tr.querySelectorAll('svg rect')].map(r => Number(r.getAttribute('width'))) })) })),
    text: s.textContent ?? '',
  }));
  const cap = INPUTS.arms.find(a => a.id === 'cap-1000').all, budget = INPUTS.arms.find(a => a.id === 'budget-1000').all;
  check(tag + ' inputs lead gives the two worse counts', read.lead.includes(cap.worse + ' of ' + cap.datasets + ' datasets') && read.lead.includes(budget.worse + ' of ' + budget.datasets), read.lead);
  const rows = read.groups.flatMap(g => g.rows);
  const full = read.head.includes('Median');
  check(tag + ' inputs table: ' + (full ? 'every column with its bars' : 'the narrow set, Change, Worse, Better and Sign p, fitting its frame'), full ? JSON.stringify(read.head) === JSON.stringify(['Change', 'Worse', 'Worse and better, drawn', 'Better', 'Of', 'Sign p', 'Median']) : JSON.stringify(read.head) === JSON.stringify(['Change', 'Worse', 'Better', 'Sign p']) && read.fits, { head: read.head, fits: read.fits });
  check(tag + ' inputs rows: every change in order with its counts, sign p' + (full ? ' and median' : ''), JSON.stringify(rows.map(r => [r.id, r.cells])) === JSON.stringify(inputsExpected.map(e => [e.id, full ? e.cells : e.narrow])), rows.map(r => [r.id, r.cells]));
  check(tag + ' inputs groups: the changed inputs, then the default run again apart', read.groups.length === 2 && read.groups[0].rows.every(r => !inputsExpected.find(e => e.id === r.id).twin) && read.groups[1].rows.length === inputsExpected.filter(e => e.twin).length && read.groups[1].rows.every(r => inputsExpected.find(e => e.id === r.id).twin), read.groups.map(g => [g.title, g.rows.map(r => r.id)]));
  // One scale for every bar: width per dataset is the same in every row.
  const scale = Math.max(...INPUTS.arms.map(a => Math.max(a.all.better, a.all.worse)));
  const perDataset = rows.flatMap(r => { const a = INPUTS.arms.find(x => x.id === r.id).all; return [[r.bars[0], a.worse], [r.bars[1], a.better]]; });
  const unit = (perDataset.find(([, n]) => n === scale) ?? [0, 1])[0] / scale;
  if (!full) check(tag + ' inputs: no bars in the narrow set', rows.every(r => r.bars.length === 0), rows.map(r => r.bars.length));
  else check(tag + ' inputs bars drawn to their counts on one scale', rows.every(r => r.bars.length === 2) && unit > 0 && perDataset.every(([w, n]) => Math.abs(w - n * unit) < 0.5), { unit, perDataset });
  check(tag + ' inputs: no internal name reaches the reader', !INTERNAL.test(read.text), read.text.match(INTERNAL)?.[0] ?? null);
  const more = section.getByRole('button', { name: 'Show each change by species' });
  await more.focus(); await page.keyboard.press('Enter');
  const opened = await section.evaluate(s => {
    const b = [...s.querySelectorAll('button')].find(x => /each change by species/.test(x.textContent));
    const region = s.querySelector('[role="region"][aria-label="Datasets better, then worse, by species"]');
    return { expanded: b?.getAttribute('aria-expanded'), controls: b?.getAttribute('aria-controls') ?? null, rows: region ? region.querySelectorAll('tbody tr').length : 0, focus: document.activeElement === b, text: s.textContent ?? '' };
  });
  check(tag + ' inputs by-species table opens by keyboard behind its control', opened.expanded === 'true' && opened.controls !== null && opened.rows === INPUTS.arms.length + 2 && opened.focus && !INTERNAL.test(opened.text), { ...opened, text: undefined });
  await page.keyboard.press('Enter');
  // The reading by learner: closed on arrival, then one group per change in
  // the main rows' order, five learners each, with the counts and shares the
  // bundled file gives; a narrow screen keeps the learner, the two counts and
  // the resistant share, and the table fits without scrolling sideways.
  const closedLearners = await section.evaluate(s => s.querySelectorAll('[data-learner-group]').length);
  const byLearner = section.getByRole('button', { name: 'Show each change by learner' });
  await byLearner.focus(); await page.keyboard.press('Enter');
  await section.locator('[data-learner-group]').first().waitFor();
  const learners = await section.evaluate(s => {
    const b = [...s.querySelectorAll('button')].find(x => /each change by learner/.test(x.textContent));
    const region = s.querySelector('[role="region"][aria-label="Each change by learner, whole lineages held out"]');
    const shown = el => getComputedStyle(el).display !== 'none';
    return {
      expanded: b?.getAttribute('aria-expanded'), controls: b?.getAttribute('aria-controls') ?? null, focus: document.activeElement === b,
      fits: region ? region.scrollWidth <= region.clientWidth : false,
      head: region ? [...region.querySelectorAll('thead th')].filter(shown).map(th => th.textContent) : [],
      groups: region ? [...region.querySelectorAll('tbody[data-learner-group]')].map(g => ({ title: g.getAttribute('data-learner-group'), rows: [...g.querySelectorAll('tr[data-learner]')].map(tr => ({ learner: tr.getAttribute('data-learner'), name: tr.querySelector('th')?.textContent ?? '', cells: [...tr.querySelectorAll('td')].filter(shown).map(td => td.textContent) })) })) : [],
      text: s.textContent ?? '',
    };
  });
  const wide = page.viewportSize().width >= 640;
  const share = v => (v * 100).toFixed(1) + '%';
  const move = xs => share(xs[0]) + ' → ' + share(xs[1]);
  const LEARNER_ORDER = ['scm', 'cart', 'rf', 'gbt', 'lasso'];
  const learnerExpected = INPUTS.arms.map(a => LEARNER_ORDER.map(l => { const r = a.learners[l]; return [l, ...(wide ? [String(r.worse), String(r.better), String(r.datasets), move(r.recall), move(r.specificity)] : [String(r.worse), String(r.better), move(r.recall)])]; }));
  const learnerSeen = learners.groups.map(g => g.rows.map(r => [r.learner, ...r.cells]));
  check(tag + ' inputs by-learner table closed on arrival, opens by keyboard behind its own control', closedLearners === 0 && learners.expanded === 'true' && learners.controls !== null && learners.focus, { closedLearners, expanded: learners.expanded, focus: learners.focus });
  check(tag + ' inputs by-learner table: 13 groups of 5 learners, each with the reading\'s counts and shares', learners.groups.length === INPUTS.arms.length && learners.groups.every(g => g.rows.length === 5) && JSON.stringify(learnerSeen) === JSON.stringify(learnerExpected), { groups: learners.groups.length, seen: learnerSeen.slice(0, 2), expected: learnerExpected.slice(0, 2) });
  check(tag + ' inputs by-learner table: ' + (wide ? 'every column' : 'the learner, the two counts and the resistant share') + ', fitting its frame on a narrow screen', learners.head.length === (wide ? 6 : 4) && (wide || learners.fits), { head: learners.head, fits: learners.fits });
  check(tag + ' inputs by-learner table: no internal name reaches the reader', !INTERNAL.test(learners.text) && !/\b(scm|cart|rf|gbt|lasso)\b/.test(learners.groups.map(g => g.rows.map(r => r.name).join(' ')).join(' ')), learners.text.match(INTERNAL)?.[0] ?? null);
  await page.keyboard.press('Enter');
}

const scrollTo = async (page, selector) => {
  await page.locator(selector).evaluate(e => e.scrollIntoView({ block: 'start' }));
  await settle(page);
};
const activeName = page => page.evaluate(() => document.activeElement?.getAttribute('aria-label') ?? document.activeElement?.textContent?.trim());
const slideNumber = async dialog => Number(((await dialog.getByRole('combobox', { name: /^Go to slide/ }).textContent()) ?? '').trim().match(/^(\d+)/)?.[1] ?? NaN);
const goSlide = async (page, dialog, index) => {
  await dialog.getByRole('button', { name: 'Close', exact: true }).focus();
  await page.keyboard.press('Home');
  for (let i = 0; i < index; i++) await page.keyboard.press('ArrowRight');
  await settle(page);
};

async function chrome(page, tag, width) {
  const c = await page.evaluate(() => {
    const root = document.getElementById('root')?.firstElementChild;
    const canvas = root?.firstElementChild;
    const bezel = getComputedStyle(root).backgroundColor, paper = getComputedStyle(canvas).backgroundColor;
    const pad = parseFloat(getComputedStyle(root).paddingTop);
    const brand = document.querySelector('header a[href="/"]');
    const mark = brand?.querySelector('svg.grm-mark');
    const word = brand?.querySelector('span');
    const nav = document.querySelector('nav[aria-label="Views"]');
    const views = nav ? [...nav.querySelectorAll('button')] : [];
    const current = views.filter(b => b.getAttribute('aria-current') === 'page');
    const rule = getComputedStyle(document.querySelector('header + main')).borderTopStyle;
    return {
      bezel, paper, pad, radius: parseFloat(getComputedStyle(canvas).borderTopLeftRadius),
      mark: mark ? mark.querySelectorAll('polygon').length : 0, motion: mark?.getAttribute('data-motion'),
      word: word?.textContent, upper: word ? getComputedStyle(word).textTransform : null,
      views: views.map(b => b.textContent), current: current.map(b => getComputedStyle(b).backgroundColor), group: nav ? getComputedStyle(nav.firstElementChild).borderTopStyle : null,
      compact: document.querySelector('nav[aria-label="View"] [role="combobox"]') !== null,
      present: (() => { const b = document.querySelector('header button[aria-label="Present"]'); return b ? { r: parseFloat(getComputedStyle(b).borderTopLeftRadius), h: b.getBoundingClientRect().height, border: getComputedStyle(b).borderTopStyle } : null; })(),
      toggle: (() => { const b = document.querySelector('header button[aria-label^="Switch to"]'); return b ? { r: parseFloat(getComputedStyle(b).borderTopLeftRadius), w: b.getBoundingClientRect().width, h: b.getBoundingClientRect().height, border: getComputedStyle(b).borderTopStyle, icons: [...b.querySelectorAll('svg')].map(s => getComputedStyle(s).opacity) } : null; })(),
      rule,
    };
  });
  check(tag + ' bezel frame and rounded canvas', c.bezel !== c.paper && c.pad >= 4 && c.pad <= 8 && c.radius >= 14, c);
  check(tag + ' brand mark and wordmark', c.mark === 7 && (width >= 640 ? c.word === 'Genomic Resistance MappingGRM' || c.word?.startsWith('Genomic Resistance Mapping') : true) && c.upper === 'uppercase', { mark: c.mark, word: c.word, upper: c.upper });
  if (c.views.length > 0) {
    check(tag + ' view nav: four views in one outlined group, one current view with a fill', JSON.stringify(c.views) === JSON.stringify(['Overview', 'Datasets', 'Species', 'Data & Methods']) && c.current.length === 1 && c.current[0] !== 'rgba(0, 0, 0, 0)' && c.group === 'solid', c);
  } else check(tag + ' compact view selector below the fit', c.compact && width < 1280, { compact: c.compact, width });
  check(tag + ' Present and the theme toggle are outlined, squared, 32px controls with one visible icon', c.present !== null && c.present.border === 'solid' && c.present.h >= 31.5 && c.present.r < c.present.h / 2 - 1 && c.toggle !== null && c.toggle.border === 'solid' && c.toggle.h >= 31.5 && c.toggle.r < c.toggle.w / 2 - 1 && c.toggle.icons.filter(o => o === '1').length === 1, { present: c.present, toggle: c.toggle });
  check(tag + ' dashed rule under the top bar', c.rule === 'dashed', c.rule);
}

/** A lone rule too long for its cell shows cut to one line and is itself
 *  the button that shows it whole, in place, without opening the row
 *  (DatasetsView's RuleCell): by click, by keyboard and, on a touch screen,
 *  by tap. Each way opens it and closes it again. Every wait is bounded, so
 *  a cell that stops working fails this check rather than the context. */
async function longRuleCheck(page, tag, ways) {
  if (longRule === undefined) { skipped.push({ check: tag + ' long rule cut to one line opens whole in place by ' + ways.join(' and '), reason: 'no lone rule of 80 characters or more in this tree' }); return; }
  await page.getByRole('searchbox', { name: 'Search datasets' }).fill(longRule.antibiotic + ' ' + longRule.species);
  const text = page.getByRole('region', { name: 'Datasets', exact: true }).getByText(longRule.scm_rules[0], { exact: true });
  const button = text.locator('xpath=ancestor::button[1]');
  const read = async () => { await settle(page); return text.evaluate(e => ({ clipped: e.scrollWidth > e.clientWidth + 1, expanded: e.closest('button')?.getAttribute('aria-expanded') ?? null, dialog: document.querySelector('dialog[open]') !== null }), null, { timeout: 2000 }).catch(() => null); };
  const act = way => (way === 'tap' ? button.tap({ timeout: 2000 }) : way === 'keyboard' ? button.focus({ timeout: 2000 }).then(() => page.keyboard.press('Enter')) : button.click({ timeout: 2000 })).then(() => true, () => false);
  await text.waitFor({ timeout: 5000 }).catch(() => {});
  // The cell measures itself once it has rendered; give the button a moment.
  await button.waitFor({ timeout: 2000 }).catch(() => {});
  const steps = [['closed', await read()]];
  for (const way of ways) for (const want of ['open', 'closed']) {
    if (steps[steps.length - 1][1]?.expanded == null) break;
    steps.push([way + ' ' + want, (await act(way)) ? await read() : null]);
  }
  const ok = steps.length === 1 + 2 * ways.length && steps.every(([name, s]) => s !== null && !s.dialog && (name.endsWith('open') ? !s.clipped && s.expanded === 'true' : s.clipped && s.expanded === 'false'));
  check(tag + ' long rule cut to one line opens whole in place by ' + ways.join(' and '), ok, { dataset: longRule.id, steps });
}

/** On a pointer device a cut-off rule names its whole text in a native
 *  tooltip (its title) and a resting pointer changes nothing: no layout, no
 *  element, no state (A34). Run after longRuleCheck, so that check's closed
 *  and keyboard reads never need the pointer moved. */
async function ruleTipCheck(page, tag) {
  const name = tag + ' cut-off rule names its whole text in a native tooltip, and hover moves nothing';
  if (longRule === undefined) { skipped.push({ check: name, reason: 'no lone rule of 80 characters or more in this tree' }); return; }
  const text = page.getByRole('region', { name: 'Datasets', exact: true }).getByText(longRule.scm_rules[0], { exact: true });
  const read = async () => { await settle(page); return text.evaluate(e => { const b = e.getBoundingClientRect(), r = e.closest('tr').getBoundingClientRect(); return { title: e.getAttribute('title'), clipped: e.scrollWidth > e.clientWidth + 1, expanded: e.closest('button')?.getAttribute('aria-expanded') ?? null, box: [b.x, b.y, b.width, b.height], row: [r.y, r.height], elements: document.getElementsByTagName('*').length, tips: document.querySelectorAll('[role="tooltip"]').length }; }, null, { timeout: 2000 }).catch(() => null); };
  // Scrolled first, so the hover itself scrolls nothing.
  await text.scrollIntoViewIfNeeded({ timeout: 2000 }).catch(() => {});
  const before = await read();
  const hovered = await text.hover({ timeout: 2000 }).then(() => true, () => false);
  const during = await read();
  await page.mouse.move(0, 0);
  const ok = before !== null && hovered && before.title === longRule.scm_rules[0] && before.clipped && before.expanded === 'false' && before.tips === 0 && JSON.stringify(before) === JSON.stringify(during);
  check(name, ok, { dataset: longRule.id, before, during });
}

/** On a touch screen a hold does what a tap does (A34): a touch held 800 ms
 *  on a cut-off rule opens it whole and a second closes it, a hold on the
 *  first of several rules opens them all, and a hold that moves 40 px does
 *  nothing. The touches go through CDP's Input.dispatchTouchEvent, as
 *  playwright's tap does, with the hold between press and lift, at a point
 *  that hits the target itself: at 390px the sticky Dataset column can lie
 *  over the rule cell, and a touch there is a touch on the row. */
async function holdCheck(page, context, tag) {
  const cdp = await context.newCDPSession(page);
  // Holds for ms, calls `during` once before the lift, and returns what it read,
  // or null when no touch reaches the element. The click Chromium sends on the
  // lift does what a tap does, so only a read before the lift shows the hold acted.
  const hold = async (locator, ms, move = 0, during = async () => null) => {
    const point = await locator.evaluate(el => {
      el.scrollIntoView({ block: 'center', inline: 'center' });
      const r = el.getBoundingClientRect(), y = r.top + r.height / 2;
      for (const x of [r.left + r.width / 2, r.left + 8, r.right - 8]) {
        const hit = x >= 0 && x < innerWidth ? document.elementFromPoint(x, y) : null;
        if (hit !== null && el.contains(hit)) return { x, y };
      }
      return null;
    }, null, { timeout: 2000 }).catch(() => null);
    if (point === null) return null;
    const { x, y } = point;
    await settle(page);
    await cdp.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x, y }] });
    if (move) { await page.waitForTimeout(50); await cdp.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: [{ x: x + move, y }] }); }
    await page.waitForTimeout(ms);
    const beforeLift = await during();
    await cdp.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
    await settle(page);
    return { beforeLift };
  };
  const region = page.getByRole('region', { name: 'Datasets', exact: true });
  const search = page.getByRole('searchbox', { name: 'Search datasets' });
  const name = tag + ' hold on a cut-off rule opens it whole, a second hold closes it, and a hold that moves does nothing';
  if (longRule === undefined) skipped.push({ check: name, reason: 'no lone rule of 80 characters or more in this tree' });
  else {
    await search.fill(longRule.antibiotic + ' ' + longRule.species);
    const text = region.getByText(longRule.scm_rules[0], { exact: true });
    await text.waitFor({ timeout: 5000 }).catch(() => {});
    const read = () => text.evaluate(e => ({ clipped: e.scrollWidth > e.clientWidth + 1, expanded: e.closest('button')?.getAttribute('aria-expanded') ?? null, dialog: document.querySelector('dialog[open]') !== null }), null, { timeout: 2000 }).catch(() => null);
    const steps = [['closed', await read(), null]];
    for (const [step, move] of [['hold open', 0], ['hold closed', 0], ['moved hold closed', 40]]) {
      const held = await hold(text, 800, move, read);
      steps.push([step, held === null ? null : await read(), held?.beforeLift ?? null]);
    }
    // Each hold's state must already hold before the lift, as well as after it.
    const as = (step, s) => s !== null && !s.dialog && (step.endsWith('open') ? !s.clipped && s.expanded === 'true' : s.clipped && s.expanded === 'false');
    const ok = steps.every(([step, after, before]) => as(step, after) && (step === 'closed' || as(step, before)));
    check(name, ok, { dataset: longRule.id, steps });
  }
  await search.fill('');
  const more = region.locator('button[aria-label^="Show all"][aria-label$="rules"]').first();
  await more.waitFor({ timeout: 5000 }).catch(() => {});
  const want = Number(((await more.getAttribute('aria-label', { timeout: 2000 }).catch(() => '')) ?? '').match(/Show all (\d+) rules/)?.[1] ?? NaN);
  // Pinned by key: once open, the row's toggle reads Fewer, and `more` moves on to the next row.
  const key = await more.evaluate(e => e.closest('tr').getAttribute('data-row-key'), null, { timeout: 2000 }).catch(() => null);
  const row = region.locator(`tr[data-row-key="${key}"]`);
  const readRow = () => row.evaluate(r => ({ rules: r.querySelectorAll('ol li').length, expanded: r.querySelector('button[aria-expanded]')?.getAttribute('aria-expanded'), dialog: document.querySelector('dialog[open]') !== null }), null, { timeout: 2000 }).catch(() => null);
  const held = await hold(row.locator('[data-row-ignore]').first(), 800, 0, readRow);
  const opened = held === null ? null : await readRow();
  const whole = s => s !== null && s.rules === want && s.expanded === 'true' && !s.dialog;
  check(tag + ' hold on the first of several rules opens them all in place', whole(opened) && whole(held?.beforeLift ?? null), { want, opened, beforeLift: held?.beforeLift ?? null });
  if (opened?.expanded === 'true') await row.getByRole('button', { name: 'Fewer' }).tap({ timeout: 2000 }).catch(() => {});
}

/** With intervals/<id>.json beside the report (`grm intervals`), the open
 *  inspector shows each random-split VME and ME the file has under its
 *  value, with the isolates it counts as a tooltip (X2). A tree without
 *  the file has no interval to show, and says so in the receipt. */
async function intervalCheck(page, tag) {
  const name = tag + ' learner table shows the offline 95% intervals of VME and ME';
  const file = path.join(FIXTURE, 'intervals', inspectId + '.json');
  if (!fs.existsSync(file)) { skipped.push({ check: name, reason: 'no intervals/' + inspectId + '.json in this tree' }); return; }
  const rows = JSON.parse(fs.readFileSync(file)).splits.find(s => s.scheme === 'random')?.learners ?? [];
  const want = rows.flatMap(r => [r.vme, r.me]).filter(v => v !== null).map(v => v.x + ' of ' + v.n);
  const shown = await page.locator('#ds-learners').evaluate(s => [...s.querySelectorAll('td span[title$="isolates called susceptible"], td span[title$="isolates called resistant"]')].map(e => ({ title: e.getAttribute('title'), text: e.textContent })), null, { timeout: 2000 }).catch(() => null);
  const counts = shown?.map(e => (e.title.match(/^([\d,]+ of [\d,]+)/)?.[1] ?? '').replace(/,/g, ''));
  const ok = shown !== null && shown.length === want.length && want.length > 0 && JSON.stringify([...counts].sort()) === JSON.stringify([...want].sort()) && shown.every(e => /^95% interval \d+\.\d%–\d+\.\d%$/.test(e.text));
  check(name, ok, { want: want.length, shown: shown?.slice(0, 4) });
}

/** The learner table must fit its frame at desktop widths with the widest
 *  interval text the tree holds, in the error view that prints it, not
 *  only for the inspected dataset's default view. */
async function widestIntervalFit(page, tag) {
  const name = tag + ' learner table with the widest interval text fits without sideways scroll';
  if (widestInterval === null) { skipped.push({ check: name, reason: 'no intervals/ in this tree' }); return; }
  const ranges = '#ds-learners td span[title$="isolates called susceptible"], #ds-learners td span[title$="isolates called resistant"]';
  await page.goto(BASE + '/datasets?dataset=' + widestInterval.id);
  await page.locator('#ds-lineages [role="listbox"]').waitFor({ timeout: 10000 }).catch(() => {});
  await page.locator(ranges).first().waitFor({ timeout: 5000 }).catch(() => {});
  if (widestInterval.scheme === 'lineage') await page.getByRole('button', { name: 'Errors, lineages held out' }).click({ timeout: 2000 }).catch(() => {});
  const frame = await page.locator('dialog [role="region"][aria-label="Learner comparison"]').evaluate(async (e, sel) => {
    await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
    return { scroll: e.scrollWidth, client: e.clientWidth, ranges: document.querySelectorAll(sel).length, pressed: document.querySelector('#ds-learners [role="group"][aria-label="Errors measured under"] [aria-pressed="true"]')?.textContent ?? null };
  }, ranges, { timeout: 5000 }).catch(() => null);
  const view = widestInterval.scheme === 'lineage' ? 'Errors, lineages held out' : 'Errors, random split';
  check(name, frame !== null && frame.ranges > 0 && frame.pressed === view && frame.scroll <= frame.client + 1, { ...widestInterval, ...frame });
}

/** At desktop widths each learner strip, axis included, ends inside its
 *  cell's content box once the table settles, in both error views, and its
 *  marks inside the cell. The cell clips a strip that is briefly too wide,
 *  so this reads the svg's own box, not what is painted. Read on the
 *  widest-interval dataset, whose lineage view narrows the strip most. */
async function stripFit(page, tag) {
  const id = widestInterval?.id ?? inspectId;
  await page.goto(BASE + '/datasets?dataset=' + id);
  await page.locator('#ds-learners table svg').first().waitFor({ timeout: 10000 }).catch(() => {});
  for (const view of ['Errors, random split', 'Errors, lineages held out']) {
    const name = tag + ' learner strips end inside their cells, ' + view;
    const button = page.locator('#ds-learners [role="group"][aria-label="Errors measured under"]').getByRole('button', { name: view, exact: true });
    if (await button.count() === 0) { skipped.push({ check: name, reason: id + ' has no "' + view + '" button' }); continue; }
    await button.click({ timeout: 2000 }).catch(() => {});
    const read = await page.locator('#ds-learners').evaluate(async (s, STATE_MS) => {
      const measure = () => {
        const cells = [...s.querySelectorAll('table tr > :nth-child(2)')].filter(c => c.querySelector('svg') !== null);
        let over = -Infinity, marksOut = 0;
        for (const cell of cells) {
          const b = cell.getBoundingClientRect(), st = getComputedStyle(cell);
          const left = b.left + parseFloat(st.paddingLeft) + parseFloat(st.borderLeftWidth), right = b.right - parseFloat(st.paddingRight) - parseFloat(st.borderRightWidth);
          for (const svg of cell.querySelectorAll('svg')) { const r = svg.getBoundingClientRect(); over = Math.max(over, r.right - right, left - r.left); }
          for (const m of cell.querySelectorAll('svg circle, svg line, svg text')) { const r = m.getBoundingClientRect(); if (r.left < b.left - 0.5 || r.right > b.right + 0.5) marksOut += 1; }
        }
        return { cells: cells.length, over: Math.round(over * 10) / 10, marksOut };
      };
      // Until the cell's observer fires and the strip redraws, up to STATE_MS.
      const t0 = performance.now();
      let m = measure();
      while (m.over > 1 && performance.now() - t0 < STATE_MS) { await new Promise(r => requestAnimationFrame(r)); m = measure(); }
      return { ...m, waited: Math.round(performance.now() - t0), pressed: s.querySelector('[role="group"][aria-label="Errors measured under"] [aria-pressed="true"]')?.textContent ?? null };
    }, STATE_MS, { timeout: STATE_MS + 5000 }).catch(() => null);
    check(name, read !== null && read.cells > 1 && read.pressed === view && read.over <= 1 && read.marksOut === 0, { id, ...read });
  }
  // Paint: every strip pushed so its right end sits 64 px past its cell's
  // content box, as a strip still drawn for a wider cell would, carries its
  // marks past the cell's edge; the Random column's pixels must not change,
  // because the cell clips them. The push is a margin, a layout change: a
  // transform painted but, in about one read in six, left the marks' boxes
  // where they were. A fixed 64 px push was not enough: at 1440 px a strip
  // ends about 46 px short of its cell's content box, so the push carried
  // its last tick about 1 px past the edge, and a pixel's change in the
  // column widths left every mark inside the cell. So the push adds each
  // cell's slack, the fonts load before the first screenshot, and the count
  // is polled each frame until every strip reads past its edge, up to
  // STATE_MS.
  const name = tag + ' a learner strip pushed past its cell paints nothing over Random';
  await page.locator('#ds-learners table').scrollIntoViewIfNeeded({ timeout: 2000 }).catch(() => {});
  await settle(page);
  const box = await page.evaluate(() => {
    const cells = [...document.querySelectorAll('#ds-learners table tr > :nth-child(3)')].map(c => c.getBoundingClientRect());
    if (cells.length === 0) return null;
    const x = Math.min(...cells.map(r => r.left)), y = Math.min(...cells.map(r => r.top));
    return { x: x + scrollX, y: y + scrollY, width: Math.max(...cells.map(r => r.right)) - x, height: Math.max(...cells.map(r => r.bottom)) - y };
  });
  if (box === null) { check(name, false, 'no Random column'); return; }
  const before = await page.screenshot({ clip: box });
  const push = await page.evaluate(async STATE_MS => {
    const svgs = [...document.querySelectorAll('#ds-learners table tr > :nth-child(2) svg')];
    const contentRight = s => { const c = s.closest('th, td'), b = c.getBoundingClientRect(), st = getComputedStyle(c); return b.right - parseFloat(st.paddingRight) - parseFloat(st.borderRightWidth); };
    const slack = svgs.map(s => Math.max(0, contentRight(s) - s.getBoundingClientRect().right));
    svgs.forEach((s, i) => { s.style.marginLeft = (slack[i] + 64) + 'px'; });
    const past = () => svgs.filter(s => { const cell = s.closest('th, td').getBoundingClientRect(); return [...s.querySelectorAll('line, circle')].some(m => m.getBoundingClientRect().right > cell.right + 1); }).length;
    const t0 = performance.now();
    let n = 0;
    do { await new Promise(r => requestAnimationFrame(r)); n = past(); } while (n < svgs.length && performance.now() - t0 < STATE_MS);
    return { moved: n, strips: svgs.length, maxSlack: Math.round(Math.max(0, ...slack) * 10) / 10, waited: Math.round(performance.now() - t0) };
  }, STATE_MS);
  const after = await page.screenshot({ clip: box });
  await page.evaluate(() => { for (const s of document.querySelectorAll('#ds-learners table tr > :nth-child(2) svg')) s.style.marginLeft = ''; });
  check(name, push.moved > 1 && before.equals(after), { id, stripsPastEdge: push.moved, strips: push.strips, maxSlack: push.maxSlack, waited: push.waited, randomUnchanged: before.equals(after) });
}

/** The widest default decision tree and the longest catalog class list the
 *  tree holds, every catalog row shown. At desktop widths each fits its
 *  frame with no sideways scroll; on a phone it may scroll inside its frame,
 *  but never the page, and its far edge scrolls into view. Every node of the
 *  tree keeps its text inside the 10 px inset its layout reserves. */
async function panelFit(page, tag, width) {
  const desktop = width >= 1280;
  for (const [key, pick, label] of [['decision tree', widestPanels.tree, 'Decision tree'], ['catalog gene groups', widestPanels.catalog, 'Catalog gene groups']]) {
    const name = tag + ' widest ' + key + (desktop ? ' fits without sideways scroll' : ' scrolls inside its frame, not the page');
    if (pick === null) { skipped.push({ check: name, reason: 'no dataset in this tree has ' + key }); continue; }
    await page.goto(BASE + '/datasets?dataset=' + pick.id);
    const region = page.locator('dialog [role="region"][aria-label="' + label + '"]');
    await region.waitFor({ timeout: 10000 }).catch(() => {});
    if (key === 'catalog gene groups') {
      const more = region.locator('xpath=ancestor::div[contains(@class, "scroll-frame")]/following-sibling::div//button');
      if (await more.count() > 0) await more.first().click({ timeout: 2000 }).catch(() => {});
    }
    await region.scrollIntoViewIfNeeded({ timeout: 2000 }).catch(() => {});
    const read = await region.evaluate(async e => {
      const frames = n => new Promise(r => { const step = k => (k === 0 ? r() : requestAnimationFrame(() => step(k - 1))); step(n); });
      await frames(2);
      const content = e.querySelector('table, svg');
      const out = { scroll: e.scrollWidth, client: e.clientWidth, page: document.documentElement.scrollWidth - innerWidth, rows: e.querySelectorAll('tbody tr').length, svg: Number(e.querySelector('svg')?.getAttribute('width') ?? 0) };
      e.scrollLeft = e.scrollWidth;
      await frames(1);
      out.reach = Math.round((content.getBoundingClientRect().right - e.getBoundingClientRect().right) * 10) / 10;
      out.labelsOut = [...e.querySelectorAll('svg g.cursor-pointer')].filter(g => { const box = g.querySelector('rect').getBoundingClientRect(); return [...g.querySelectorAll('text')].some(t => { const r = t.getBoundingClientRect(); return r.left < box.left - 0.5 || r.right > box.right - 10 + 0.5; }); }).length;
      e.scrollLeft = 0;
      return out;
    }, null, { timeout: 5000 }).catch(() => null);
    // The pick is the one on screen: the tree at its column count (at least
    // 140 px each), every catalog row.
    const shown = read !== null && (key === 'decision tree' ? read.svg >= pick.columns * 140 + 8 : read.rows === pick.groups);
    const fits = read !== null && (desktop ? read.scroll <= read.client + 1 : read.page <= 1 && read.reach <= 1);
    check(name, shown && fits, { ...pick, ...read });
    if (key !== 'decision tree') continue;
    // Read again once the fonts settle and the tree has measured again: every
    // label still inside its node, and no leaf ("Calls resistant", "Calls
    // susceptible", which fit any node in Inter) cut short.
    const settled = await region.evaluate(async e => {
      await document.fonts.ready;
      await new Promise(r => setTimeout(r, 250));
      await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
      const nodes = [...e.querySelectorAll('svg g.cursor-pointer')];
      const labelsOut = nodes.filter(g => { const box = g.querySelector('rect').getBoundingClientRect(); return [...g.querySelectorAll('text')].some(t => { const r = t.getBoundingClientRect(); return r.left < box.left - 0.5 || r.right > box.right - 10 + 0.5; }); }).length;
      const leavesCut = nodes.filter(g => /^Calls .*…$/.test(g.querySelector('text')?.textContent ?? '')).length;
      return { labelsOut, leavesCut, fonts: document.fonts.status };
    }, null, { timeout: 5000 }).catch(() => null);
    check(tag + ' widest decision tree keeps each label inside its node\'s text inset', read !== null && read.labelsOut === 0 && settled !== null && settled.labelsOut === 0 && settled.leavesCut === 0, { id: pick.id, labelsOut: read?.labelsOut, settled });
  }
  await labelFit(page, tag);
}

/** The decision tree's labels in widestPanels.labels, once the fonts settle:
 *  every line ends inside its node's 10 px text inset, and a line is cut
 *  (…) or a count line shortened (no "isolates") only when its full text,
 *  measured in a copy of the drawn line, is wider than that room. A split's
 *  full label is read from the tree's own listbox option. */
async function labelFit(page, tag) {
  const name = tag + ' decision tree labels fit their nodes, cut only when the full text does not fit';
  if (widestPanels.labels.length === 0) { skipped.push({ check: name, reason: 'no dataset in this tree has a decision tree' }); return; }
  const trees = [];
  for (const id of widestPanels.labels) {
    await page.goto(BASE + '/datasets?dataset=' + id);
    const region = page.locator('dialog [role="region"][aria-label="Decision tree"]');
    await region.waitFor({ timeout: 10000 }).catch(() => {});
    await region.scrollIntoViewIfNeeded({ timeout: 2000 }).catch(() => {});
    await region.evaluate(async () => { await document.fonts.ready; await new Promise(r => setTimeout(r, 250)); await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r))); }, null, { timeout: 5000 }).catch(() => {});
    // Every node's option, by the keyboard path, in the nodes' order.
    const listbox = region.locator('[role="listbox"]');
    const nodes = await region.evaluate(e => e.querySelectorAll('svg g.cursor-pointer').length).catch(() => 0);
    const options = [];
    if (await listbox.count() > 0) {
      await listbox.focus().catch(() => {});
      for (let i = 0; i < nodes; i++) {
        await page.keyboard.press(i === 0 ? 'Home' : 'ArrowRight');
        options.push(await listbox.evaluate(e => e.querySelector('[role="option"]')?.textContent ?? null).catch(() => null));
      }
      await page.keyboard.press('Escape');
      await listbox.evaluate(e => e.blur()).catch(() => {});
    }
    const read = await region.evaluate((e, options) => {
      const bad = [];
      let lines = 0, cut = 0;
      [...e.querySelectorAll('svg g.cursor-pointer')].forEach((g, i) => {
        const box = g.querySelector('rect').getBoundingClientRect(), room = box.width - 20;
        [...g.querySelectorAll('text')].forEach((t, j) => {
          lines += 1;
          const r = t.getBoundingClientRect(), text = t.textContent;
          if (r.left < box.left - 0.5 || r.right > box.right - 10 + 0.5) bad.push({ node: i, text, past: Math.round((r.right - box.right + 10) * 10) / 10 });
          const short = j === 1 && !text.includes(' isolates · ');
          if (!text.endsWith('…') && !short) return;
          cut += 1;
          let full = null;
          if (j === 1) full = text.replace(' · ', ' isolates · ');
          else if (/^Calls r/.test(text)) full = 'Calls resistant';
          else if (/^Calls s/.test(text)) full = 'Calls susceptible';
          else full = /^Split on (.*), [\d,]+ training isolates, /.exec(options[i] ?? '')?.[1] ?? null;
          if (full === null) { bad.push({ node: i, text, full: 'unknown' }); return; }
          const copy = t.cloneNode(false);
          copy.textContent = full;
          t.parentNode.appendChild(copy);
          const w = copy.getComputedTextLength();
          copy.remove();
          if (w + 0.5 <= room) bad.push({ node: i, text, full, fullWidth: Math.round(w * 10) / 10, room });
        });
      });
      return { lines, cut, bad };
    }, options, { timeout: 5000 }).catch(() => null);
    trees.push({ id, nodes, options: options.filter(o => o !== null).length, ...(read ?? { lines: 0, cut: 0, bad: [{ error: 'unread' }] }) });
  }
  const ok = trees.every(t => t.lines > 0 && t.options === t.nodes && t.bad.length === 0);
  check(name, ok, trees.map(t => ({ id: t.id, lines: t.lines, cut: t.cut, options: t.options, nodes: t.nodes, bad: t.bad.length, first: t.bad.slice(0, 3) })));
}

async function interactions(page, tag, width) {
  await page.goto(BASE + '/datasets');
  const search = page.getByRole('searchbox', { name: 'Search datasets' });
  await search.waitFor();
  await chrome(page, tag, width);
  await page.locator('header a[href="/"]').focus();
  const focus = [];
  for (let i = 0; i < 7; i++) { focus.push(await activeName(page)); await page.keyboard.press('Tab'); }
  const brandName = width >= 640 ? 'Genomic Resistance MappingGRM' : 'Genomic Resistance MappingGRM';
  const theme = tag.includes('dark') ? 'Switch to light theme' : 'Switch to dark theme';
  const wanted = width >= 1280 ? [brandName, 'Overview', 'Datasets', 'Species', 'Data & Methods', 'Present', theme] : [brandName, 'View: Datasets', 'Present', theme];
  check(tag + ' top bar keyboard order', JSON.stringify(focus.slice(0, wanted.length)) === JSON.stringify(wanted), focus);
  if (width >= 1280) check(tag + ' current view marked', await page.locator('nav[aria-label="Views"] [aria-current="page"]').textContent() === 'Datasets', null);
  await headerCheck(page, tag + ' datasets', 'Datasets');
  await controlCheck(page, tag + ' datasets');
  await edgeCheck(page, tag + ' datasets', '[data-view-scroll] h1');
  const scopeText = await page.locator('[data-view-scroll] header').textContent();
  check(tag + ' plain scope sentence', scopeText.includes(CHIP) && !scopeText.includes('BV-BRC lab records'), scopeText);
  await shot(page, tag + '-status-line');
  await page.mouse.move(0, 0);
  const inline = width >= 1440;
  if (inline) {
    const align = await page.evaluate(() => {
      const aside = document.querySelector('aside[data-view-options-mode="inline"]');
      const row = aside.firstElementChild.getBoundingClientRect(), input = document.querySelector('input[type="search"]').getBoundingClientRect();
      const group = aside.querySelector('.filter-group-label').getBoundingClientRect(), th = document.querySelector('[role="region"][aria-label="Datasets"] thead').getBoundingClientRect();
      return { rowTop: row.top, rowHeight: row.height, inputTop: input.top, inputHeight: input.height, groupTop: group.top, headTop: th.top };
    });
    check(tag + ' rail title row on the search row', Math.abs(align.rowTop - align.inputTop) <= 1 && Math.abs(align.rowHeight - align.inputHeight) <= 1, align);
    check(tag + ' rail first group level with the table header', Math.abs(align.groupTop - align.headTop) <= 8, align);
    const toggle = page.getByRole('button', { name: 'Collapse filters' });
    const before = await toggle.evaluate(e => e.getBoundingClientRect().toJSON());
    await toggle.click(); await settle(page);
    const collapsed = await page.evaluate(() => { const a = document.querySelector('aside[data-view-options-mode="inline"]'); const b = a.querySelector('button'); return { width: a.getBoundingClientRect().width, toggle: b.getBoundingClientRect().toJSON(), name: b.getAttribute('aria-label') }; });
    check(tag + ' collapsed strip 3.5rem, toggle on the same pixel', Math.abs(collapsed.width - 57) <= 1 && Math.abs(collapsed.toggle.x - before.x) <= 0.5 && Math.abs(collapsed.toggle.y - before.y) <= 0.5, { before, collapsed });
    await shot(page, tag + '-rail-collapsed');
    await page.getByRole('button', { name: /^Expand filters/ }).click(); await settle(page);
  } else {
    await page.getByRole('button', { name: 'Filters', exact: true }).click();
  }
  await page.getByRole('checkbox', { name: /^Has MIC prediction/ }).click();
  await shot(page, tag + '-filters');
  const band = await page.locator('[data-view-options-panel] footer').textContent();
  check(tag + ' rail band: the filter count and Clear once a filter is set', band === '1 filterClear', band);
  const bandEdge = await page.evaluate(() => { const panel = document.querySelector('[data-view-options-panel]'); return { clear: panel.querySelector('footer button')?.getBoundingClientRect().right, field: panel.querySelector('.filter-field-control [role="combobox"]')?.getBoundingClientRect().right }; });
  check(tag + ' rail band: Clear ends on the right edge of the controls', bandEdge.clear !== undefined && bandEdge.field !== undefined && Math.abs(bandEdge.clear - bandEdge.field) <= 1, bandEdge);
  const groupGaps = await page.evaluate(() => [...document.querySelectorAll('[data-view-options-panel] .filter-rail-body > * + *')].map(e => getComputedStyle(e).marginTop));
  check(tag + ' rail groups share one spacing', groupGaps.length > 0 && new Set(groupGaps).size === 1, groupGaps);
  if (!inline) {
    await page.keyboard.press('Escape');
    await page.locator('aside[data-view-options-mode="sheet"]').waitFor({ state: 'detached' });
    await page.waitForFunction(() => document.activeElement?.hasAttribute('data-filter-toggle'));
    check(tag + ' filter focus return', (await activeName(page)).startsWith('Filters'), await activeName(page));
  }
  check(tag + ' MIC filter count in the search row', await page.getByText(expected.mic.length + ' of ' + expected.datasets + ' datasets', { exact: true }).count() >= 1, null);
  if (inline) await shot(page, tag + '-filters-docked-active');
  if (!inline) await page.getByRole('button', { name: /^Filters/ }).click();
  await page.locator('[data-view-options-panel] footer').getByRole('button', { name: 'Clear', exact: true }).click();
  if (!inline) {
    await page.keyboard.press('Escape');
    await page.locator('aside[data-view-options-mode="sheet"]').waitFor({ state: 'detached' });
  }
  check(tag + ' filters cleared', await page.getByText(expected.datasets + ' datasets', { exact: true }).count() >= 1, null);
  await search.fill('zz-no-dataset-matches-zz');
  await page.getByRole('button', { name: 'Clear filters and search' }).waitFor();
  await shot(page, tag + '-empty-filter');
  await page.getByRole('button', { name: 'Clear filters and search' }).click();
  check(tag + ' clear filter recovery', await search.inputValue() === '', null);
  // The table offers Show all only past one page of rows (PAGE_SIZE in DataTable.tsx).
  const showAll = page.getByRole('button', { name: 'Show all ' + expected.datasets.toLocaleString('en-US'), exact: true });
  if (expected.datasets > 25) await showAll.click();
  else check(tag + ' one page of rows offers no Show all', await showAll.count() === 0, null);
  const status = await page.getByRole('region', { name: 'Datasets', exact: true }).evaluate(region => {
    const rows = [...region.querySelectorAll('tbody tr')];
    const v = (r, i) => r.children[i]?.querySelector('[data-verdict]')?.getAttribute('data-verdict');
    return { rows: rows.length, vme: rows.filter(r => v(r, 5) === 'meets').length, me: rows.filter(r => v(r, 6) === 'meets').length, mic: rows.filter(r => v(r, 8) === 'meets').length, headers: [...region.querySelectorAll('thead th')].map(e => e.textContent), scope: [...region.querySelectorAll('thead th')].every(e => e.scope === 'col') };
  });
  check(tag + ' judged thresholds and counts', status.rows === expected.datasets && status.vme === expected.stats.vmeMeets && status.me === expected.stats.meMeets && status.mic === expected.stats.micAbove90 && status.scope, status);
  check(tag + ' table headers', JSON.stringify(status.headers) === JSON.stringify(['Dataset', 'Random', 'Lineage', 'Gappoints', 'Isolates', 'VME≤1.5%', 'ME≤3%', 'Catalog', 'MIC≥90%', 'Ruleset covering machine']), status.headers);
  await judgedCheck(page, tag);
  const classes = await page.getByRole('region', { name: 'Datasets', exact: true }).evaluate(region => {
    const rows = [...region.querySelectorAll('tbody tr')];
    return {
      rows: rows.length,
      marks: rows.map(r => { const sub = r.querySelector('.dataset-sub'); const g = sub?.querySelector('svg.class-glyph'); const shape = g?.querySelector('circle, rect, path'); const name = (sub?.textContent ?? '').split('· ').pop(); return { name, family: g?.getAttribute('data-class-family') ?? null, shape: g?.getAttribute('data-shape') ?? null, hidden: g?.getAttribute('aria-hidden'), fill: shape ? getComputedStyle(shape).fill : null, stroke: shape ? getComputedStyle(shape).stroke : null }; }),
      key: document.querySelectorAll('[data-class-key]').length,
    };
  });
  const byFamily = new Map();
  for (const m of classes.marks) { const paint = m.shape === 'ring' ? m.stroke : m.fill; if (!byFamily.has(m.family)) byFamily.set(m.family, new Set()); byFamily.get(m.family).add(m.shape + ' ' + paint); }
  const paints = [...byFamily.values()].map(v => [...v]);
  check(tag + ' every row names its drug class beside a decorative family mark', classes.rows === expected.datasets && classes.marks.every(m => m.family !== null && m.hidden === 'true' && m.name && manifest.datasets.some(d => d.drug_class === m.name)), classes.marks.filter(m => m.family === null || m.hidden !== 'true'));
  check(tag + ' one shape and one paint per family, distinct across families, with one key', classes.key === 1 && paints.every(p => p.length === 1) && new Set(paints.map(p => p[0])).size === paints.length, { key: classes.key, paints: Object.fromEntries([...byFamily.entries()].map(([k, v]) => [k, [...v]])) });
  // Off every header first: a hovered header is deep blue, as it should be.
  await page.mouse.move(2, 2);
  const heads = await page.getByRole('region', { name: 'Datasets', exact: true }).evaluate(region => {
    const blue = getComputedStyle(document.documentElement).getPropertyValue('--blue').trim().split(/\s+/).join(', ');
    return [...region.querySelectorAll('thead th')].map(th => { const label = th.querySelector('[data-th-label]'); const glyph = th.querySelector('svg'); return { text: th.textContent, color: getComputedStyle(label).color, blue: 'rgb(' + blue + ')', glyph: glyph ? Number(getComputedStyle(glyph).opacity) : null, sortable: th.querySelector('button') !== null }; });
  });
  check(tag + ' blue header words', heads.every(h => h.color === h.blue), heads.map(h => [h.text, h.color, h.blue]));
  check(tag + ' sort glyph slot always occupied', heads.filter(h => h.sortable).every(h => h.glyph !== null && h.glyph >= 0.4), heads.map(h => [h.text, h.glyph]));
  const align = await page.getByRole('region', { name: 'Datasets', exact: true }).evaluate(region => {
    const th = [...region.querySelectorAll('thead th')];
    const row = region.querySelector('tbody tr');
    return [1, 2, 4, 7].map(i => { const label = th[i].querySelector('[data-th-label]'); const cell = row.children[i]; const pad = parseFloat(getComputedStyle(cell).paddingRight); return { header: th[i].textContent, labelRight: label.getBoundingClientRect().right, valueRight: cell.getBoundingClientRect().right - pad }; });
  });
  check(tag + ' numeric headers line up with values', align.every(a => Math.abs(a.labelRight - a.valueRight) <= 1.5), align);
  const table = page.getByRole('region', { name: 'Datasets', exact: true });
  // Species and class metadata may use a second line; a phone may also
  // wrap a long antibiotic name. Keep both bounded without requiring the
  // desktop and phone to have identical row heights.
  const rowMeta = await table.evaluate(region => [...region.querySelectorAll('tbody tr')].slice(0, 25).map(r => {
    const sub = r.querySelector('.dataset-sub');
    const line = sub ? parseFloat(getComputedStyle(sub).lineHeight) : 0;
    return { height: Math.round(r.getBoundingClientRect().height), subtitleLines: sub && line > 0 ? sub.getBoundingClientRect().height / line : 99 };
  }));
  const heights = rowMeta.map(r => r.height);
  const shortest = Math.min(...heights), tallest = Math.max(...heights);
  const maxHeight = width < 640 ? 100 : 80, maxSpread = width < 640 ? 40 : 20;
  check(tag + ' compact rows', shortest >= 44 && tallest <= maxHeight && tallest - shortest <= maxSpread && rowMeta.every(r => r.subtitleLines <= 2.1), rowMeta);
  const more = table.locator('button[aria-label^="Show all"][aria-label$="rules"]').first();
  // Pin the row itself: once expanded, the first "+N" is the next row's.
  const ruleRow = await more.evaluateHandle(e => e.closest('tr'));
  // Its row has had focus (the button), so the unfocused-row read skips it.
  const ruleKey = await ruleRow.evaluate(r => r.getAttribute('data-row-key'));
  const shortH = await ruleRow.evaluate(r => r.getBoundingClientRect().height);
  const wantRules = Number((await more.getAttribute('aria-label')).match(/Show all (\d+) rules/)[1]);
  await more.click(); await settle(page);
  const opened = await ruleRow.evaluate(r => ({ height: r.getBoundingClientRect().height, rules: r.querySelectorAll('ol li').length, expanded: r.querySelector('button[aria-expanded="true"]')?.textContent }));
  check(tag + ' rule expands in place', opened.rules === wantRules && opened.height > shortH && opened.expanded === 'Fewer' && !(await page.locator('dialog[open]').count()), { shortH, opened, wantRules });
  await shot(page, tag + '-rule-open');
  await table.getByRole('button', { name: 'Fewer' }).click();
  const left = await table.locator('tbody tr').first().locator('td').first().evaluate(e => e.getBoundingClientRect().left);
  await table.evaluate(e => e.scrollLeft = 500);
  const shifted = await table.locator('tbody tr').first().locator('td').first().evaluate(e => e.getBoundingClientRect().left);
  check(tag + ' sticky dataset column', Math.abs(left - shifted) <= 1.5, { left, shifted });
  await shot(page, tag + '-dataset-status');
  await table.evaluate(e => e.scrollLeft = 0);
  await page.getByRole('button', { name: 'Random', exact: true }).focus();
  await page.keyboard.press('Enter');
  check(tag + ' keyboard sorting', await table.locator('th[aria-sort="descending"]').textContent() === 'Random', null);
  const rows = table.locator('[data-row-open]');
  // Off the table first: a hovered row has a tint of its own (blue at 3.5%),
  // which must not stand in for the focus tint.
  await page.mouse.move(2, 2);
  await rows.nth(0).focus(); await page.keyboard.press('ArrowDown');
  check(tag + ' keyboard row movement', await rows.nth(1).evaluate(e => e === document.activeElement), null);
  // The tint is a colour transition (transition-colors on the row), so it
  // lands a frame or more after focus moves. Read it after two frames, then
  // each frame until two frames in a row agree on a colour that isn't
  // clear, for up to STATE_MS, and compare that colour with the focus tint,
  // blue at 5% (DataTable's focus-within:bg-blue/[0.05]).
  const tint = await rows.nth(1).evaluate((e, STATE_MS) => new Promise(resolve => {
    const row = e.closest('tr'), start = performance.now();
    const probe = document.createElement('div');
    probe.style.backgroundColor = 'rgb(' + getComputedStyle(document.documentElement).getPropertyValue('--blue').trim() + ' / 0.05)';
    document.body.append(probe); const want = getComputedStyle(probe).backgroundColor; probe.remove();
    let last = null;
    const read = () => { const c = getComputedStyle(row).backgroundColor; const settled = c === last && c !== 'rgba(0, 0, 0, 0)'; if (settled || performance.now() - start > STATE_MS) resolve({ colour: c, settled, want, ms: Math.round(performance.now() - start) }); else { last = c; requestAnimationFrame(read); } };
    requestAnimationFrame(() => requestAnimationFrame(read));
  }), STATE_MS);
  check(tag + ' focused row tints', tint.settled && tint.colour === tint.want, tint);
  // The tint follows focus, in the row and in its cells. The sticky cell is
  // left out: it paints its own paper (and its tint has a check of its
  // own). A colour with alpha 0 is clear. With the pointer still off the
  // table, a row that never had focus (not the first two rows, which
  // ArrowDown used, nor the row whose rules were opened) has no tint...
  const unfocused = await table.evaluate((region, ruleKey) => {
    const clear = c => c === 'transparent' || /^rgba\([^,]+,[^,]+,[^,]+,\s*0\)$/.test(c);
    const rs = [...region.querySelectorAll('tbody tr')];
    const index = rs.findIndex((r, k) => k >= 2 && r.getAttribute('data-row-key') !== ruleKey);
    const r = rs[index];
    if (r === undefined) return { index, tinted: null };
    const paints = [['tr', getComputedStyle(r).backgroundColor], ...[...r.children].map((td, k) => ['td' + k, td]).filter(([, td]) => getComputedStyle(td).position !== 'sticky').map(([k, td]) => [k, getComputedStyle(td).backgroundColor])];
    return { index, key: r.getAttribute('data-row-key'), tinted: paints.filter(([, c]) => !clear(c)) };
  }, ruleKey);
  check(tag + ' unfocused row has no tint', unfocused.tinted !== null && unfocused.tinted.length === 0, unfocused);
  // ... and the row ArrowDown left loses its tint, read each frame for up to
  // STATE_MS until it is clear.
  const leftRow = await rows.nth(0).evaluate((e, STATE_MS) => new Promise(resolve => {
    const r = e.closest('tr'), start = performance.now();
    const clear = c => c === 'transparent' || /^rgba\([^,]+,[^,]+,[^,]+,\s*0\)$/.test(c);
    const read = () => {
      const paints = [['tr', getComputedStyle(r).backgroundColor], ...[...r.children].map((td, k) => ['td' + k, td]).filter(([, td]) => getComputedStyle(td).position !== 'sticky').map(([k, td]) => [k, getComputedStyle(td).backgroundColor])];
      const tinted = paints.filter(([, c]) => !clear(c));
      if (tinted.length === 0 || performance.now() - start > STATE_MS) resolve({ tinted, ms: Math.round(performance.now() - start) }); else requestAnimationFrame(read);
    };
    requestAnimationFrame(read);
  }), STATE_MS);
  check(tag + ' tint leaves the row focus left', leftRow.tinted.length === 0, leftRow);
  const stickyTint = await table.evaluate(region => { const rs = [...region.querySelectorAll('tbody tr')]; return { focused: getComputedStyle(rs[1].children[0]).backgroundColor, rest: getComputedStyle(rs[3].children[0]).backgroundColor }; });
  check(tag + ' focus tint reaches the sticky cell', stickyTint.focused !== stickyTint.rest, stickyTint);
  await shot(page, tag + '-row-focus');
  // The app keeps a launch point for FRESH_MS (launch.ts) and the inspector
  // reads it once, in the layout effect of its first render. So the check
  // times the open itself: from the Enter keydown (a window listener in the
  // capture phase, which runs before the app's own document listener records
  // the point) to the panel's insertion (a mutation observer, called after
  // the commit that ran the layout effect). That span is never shorter than
  // the app's own. An open without a launch point is a failure, unless the
  // span reached FRESH_MS: then the machine, not the app, missed the window,
  // and the check closes the inspector and opens it again, up to three times.
  // A launch point missed inside the window fails at once.
  const openFromRow = async () => {
    await page.evaluate(() => {
      const t = window.__grmOpen = { key: null, open: null };
      window.addEventListener('keydown', e => { if (e.key === 'Enter' && t.key === null) t.key = performance.now(); }, { capture: true, once: true });
      const seen = () => document.querySelector('dialog[open] .inspector-in') !== null;
      const mo = new MutationObserver(() => { if (seen()) { t.open = performance.now(); mo.disconnect(); } });
      mo.observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ['open'] });
    });
    await page.keyboard.press('Enter');
    await page.locator('#ds-learners').waitFor();
    // Focus moves to the inspector's body in an effect after the open; wait
    // for it rather than read it on the first frame.
    await page.waitForFunction(() => document.activeElement === document.querySelector('dialog[open] [data-modal-scroll-owner]'), null, { timeout: STATE_MS }).catch(() => {});
    return page.evaluate(() => {
      const panel = document.querySelector('dialog[open] .inspector-in');
      const owner = document.querySelector('dialog[open] [data-modal-scroll-owner]');
      const t = window.__grmOpen;
      return { origin: panel ? getComputedStyle(panel).transformOrigin : null, launchY: panel?.style.getPropertyValue('--launch-y'), focus: document.activeElement === owner, tips: document.querySelectorAll('[role="tooltip"]').length, openMs: t.key !== null && t.open !== null ? Math.round(t.open - t.key) : null };
    });
  };
  const openings = [];
  for (;;) {
    const o = await openFromRow();
    openings.push(o);
    const missedWindow = o.launchY === '' && o.openMs !== null && o.openMs >= FRESH_MS;
    if (!missedWindow || openings.length === 3) break;
    await page.keyboard.press('Escape');
    await page.waitForFunction(() => !document.querySelector('dialog[open]'));
    await rows.nth(1).focus();
  }
  const opening = openings.at(-1);
  check(tag + ' inspector opens from the row, focus on its body, no tip', opening.launchY !== '' && opening.focus && opening.tips === 0, { ...opening, freshMs: FRESH_MS, attempts: openings.map(o => ({ launchY: o.launchY, openMs: o.openMs })) });
  await page.keyboard.press('Escape');
  await page.waitForFunction(() => !document.querySelector('dialog[open]'));
  check(tag + ' inspector focus return', await rows.nth(1).evaluate(e => e === document.activeElement), await activeName(page));
  await longRuleCheck(page, tag, ['click', 'keyboard']);
  await ruleTipCheck(page, tag);

  await page.goto(BASE + '/datasets?dataset=' + inspectId);
  await page.locator('#ds-model').waitFor();
  // The map is a separate, larger request. Wait for it before comparing
  // section positions and scroll-spy labels on slower public connections.
  await page.locator('#ds-lineages [role="listbox"]').waitFor();
  const head = await page.evaluate(() => {
    const header = document.querySelector('dialog[open] header');
    return { status: header.querySelector('p')?.textContent, links: [...header.querySelectorAll('nav[aria-label="Sections"] button')].map(b => b.textContent), current: header.querySelector('nav[aria-label="Sections"] [aria-current="true"]')?.textContent, title: document.getElementById('inspector-title')?.getBoundingClientRect().x, body: document.querySelector('#ds-learners h2')?.getBoundingClientRect().x, linksX: header.querySelector('nav[aria-label="Sections"] button')?.getBoundingClientRect().x, classMark: header.querySelector('p svg.class-glyph')?.getAttribute('data-class-family') ?? null };
  });
  check(tag + ' inspector header: chip and facts in one status row', head.status?.includes(CHIP.replace('Provisional metrics', 'Provisional')) || head.status?.includes('Tested within one public database only') || head.status?.includes('Provisional') || head.status?.includes('Unaudited, provisional'), head.status);
  check(tag + ' inspector section links', JSON.stringify(head.links) === JSON.stringify(['Learners', 'Interpretable model', 'Genes and panels', 'Marker designer', 'MIC', 'Lineages']) && head.current === 'Learners', head);
  check(tag + ' inspector header and body share a column', Math.abs(head.title - head.body) <= 1 && Math.abs(head.linksX - head.body) <= 1, head);
  check(tag + ' inspector names the drug class with its family mark', head.classMark !== null, head.classMark);
  const seg = await page.getByRole('group', { name: 'Errors measured under' }).evaluate(g => [...g.querySelectorAll('button')].map(b => ({ text: b.textContent, pressed: b.getAttribute('aria-pressed'), fill: getComputedStyle(b).backgroundColor, ink: 'rgb(' + getComputedStyle(document.documentElement).getPropertyValue('--ink').trim().split(/\s+/).join(', ') + ')' })));
  check(tag + ' ink segmented errors switch', seg.length === 2 && seg[0].text === 'Errors, random split' && seg[1].text === 'Errors, lineages held out' && seg[0].pressed === 'true' && seg[0].fill === seg[0].ink, seg);
  if (width >= 1280) {
    const frames = await page.evaluate(() => [...document.querySelectorAll('dialog [role="region"]')].filter(e => e.getClientRects().length > 0).map(e => ({ label: e.getAttribute('aria-label'), scroll: e.scrollWidth, client: e.clientWidth })).filter(f => f.scroll > f.client + 1));
    check(tag + ' inspector tables fit without sideways scroll', frames.length === 0, frames);
  }
  for (const id of ['ds-model', 'ds-mechanisms', 'ds-designer', 'ds-mic', 'ds-lineages']) {
    const want = { 'ds-model': 'Interpretable model', 'ds-mechanisms': 'Genes and panels', 'ds-designer': 'Marker designer', 'ds-mic': 'MIC', 'ds-lineages': 'Lineages' }[id];
    // The link follows the scroll through an observer; wait for it, up to
    // STATE_MS, then read whatever link is current.
    await scrollTo(page, '#' + id);
    await page.waitForFunction(w => document.querySelector('dialog[open] nav[aria-label="Sections"] [aria-current="true"]')?.textContent === w, want, { timeout: STATE_MS }).catch(() => {});
    const current = await page.locator('dialog[open] nav[aria-label="Sections"] [aria-current="true"]').textContent();
    check(tag + ' section link follows ' + id, current === want, current);
    await shot(page, tag + '-' + id);
  }
  await page.getByRole('button', { name: 'Method and coverage', exact: true }).click();
  await reachCheck(page, tag + ' inspector');
  await page.getByText('Evaluation audit', { exact: true }).last().scrollIntoViewIfNeeded();
  await shot(page, tag + '-method-fold');
  const fold = await page.evaluate(() => { const dt = [...document.querySelectorAll('dialog dt')].find(e => e.textContent === 'Evaluation audit'); const dd = dt?.nextElementSibling; return { text: dd?.textContent, align: dd ? getComputedStyle(dd).textAlign : null }; });
  check(tag + ' audit row reads left-aligned without undefined', fold.text !== undefined && !fold.text.includes('undefined') && fold.align !== 'right', fold);
  const owners = await page.getByRole('dialog').evaluate(d => [...d.querySelectorAll('*')].filter(e => getComputedStyle(e).overflowY === 'auto' && e.scrollHeight > e.clientHeight + 1).map(e => e.getAttribute('data-modal-scroll-owner')));
  check(tag + ' inspector vertical scroll owner', owners.length === 1 && owners[0] === 'true', owners);
  await page.keyboard.press('Escape');

  await page.goto(BASE + '/species/' + speciesSlug);
  await page.getByRole('grid', { name: /lineage block/ }).first().waitFor();
  await headerCheck(page, tag + ' species', speciesName);
  await edgeCheck(page, tag + ' species', '[data-view-scroll] h1');
  const purpose = await page.evaluate(() => document.querySelector('[data-view-scroll] header h1 + p')?.textContent);
  check(tag + ' species purpose line', purpose === "A lineage is a group of closely related genomes, like one branch of a family tree. This page shows how one species' genomes fall into lineages, and how resistance follows them.", purpose);
  const drug = page.getByRole('combobox', { name: /^Map drug/ });
  // The outline must sit on its drug's own row. The map's label column
  // draws each drug's name in a row of its own (Structure.tsx), so the
  // outline's middle must be within 2 px of the chosen drug's label's
  // middle, before the choice and after it; no height is assumed. Between the
  // two reads it must move as far as the label does, within 0.5 px, so an
  // outline that steps rows at another pitch than the labels fails. The choice
  // reaches the outline through a state update and a transform transition,
  // so read after two frames, then each frame until two frames in a row
  // agree on a position on that row, for up to STATE_MS. A fast machine
  // settles in a few frames, as before; the deadline only bounds a slow one.
  const outlineOn = label => page.evaluate(([label, STATE_MS]) => new Promise(resolve => {
    const place = () => {
      const o = [...document.querySelectorAll('[role="grid"] ~ div[aria-hidden="true"]')].find(e => e.className.includes('border-ink'));
      const l = o && [...(o.closest('#sp-structure') ?? document).querySelectorAll('span.whitespace-nowrap')].find(s => s.textContent === label && s.className.includes('leading-[16px]'));
      if (!o || !l) return null;
      const a = o.getBoundingClientRect(), b = l.getBoundingClientRect();
      return { off: Math.round((a.top + a.bottom - b.top - b.bottom) * 50) / 100, outline: [a.top, a.height], labelBox: [b.top, b.height] };
    };
    const start = performance.now();
    let last = null;
    const poll = () => {
      const p = place(), settled = p !== null && last !== null && Math.abs(p.outline[0] - last.outline[0]) < 0.01;
      if ((settled && Math.abs(p.off) <= 2) || performance.now() - start > STATE_MS) resolve({ label, ...p, settled, ms: Math.round(performance.now() - start) }); else { last = p; requestAnimationFrame(poll); }
    };
    requestAnimationFrame(() => requestAnimationFrame(poll));
  }), [label, STATE_MS]);
  await drug.click();
  const options = await page.getByRole('option').evaluateAll(os => os.map(o => ({ label: (o.querySelector('span.flex-1') ?? o).textContent.trim(), selected: o.getAttribute('aria-selected') === 'true' })));
  const fromLabel = options.find(o => o.selected)?.label ?? null, toLabel = options[1]?.label ?? null;
  const was = await outlineOn(fromLabel);
  await page.getByRole('option').nth(1).click();
  const now = await outlineOn(toLabel);
  const onRow = p => p.settled && Math.abs(p.off) <= 2;
  check(tag + ' chosen drug outline moves to its row', fromLabel !== null && toLabel !== null && fromLabel !== toLabel && onRow(was) && onRow(now)
    && Math.abs((now.outline[0] - was.outline[0]) - (now.labelBox[0] - was.labelBox[0])) <= 0.5, { was, now });
  await scrollTo(page, '#sp-cross');
  await shot(page, tag + '-cross-matrix');
  await shot(page, tag + '-cross-table');
  check(tag + ' cross-resistance table', await page.getByRole('region', { name: 'Cross-resistance pairs' }).count() === 1, null);
  await page.goto(BASE + '/tests');
  await page.getByRole('region', { name: 'Panel results' }).waitFor();
  check(tag + ' legacy tests redirect', new URL(page.url()).pathname === '/' && new URL(page.url()).hash === '#panels', page.url());
  await shot(page, tag + '-panel-evidence');
  await page.goto(BASE + '/methods');
  await page.getByRole('heading', { name: /Models and evaluation/ }).waitFor();
  await headerCheck(page, tag + ' methods', 'Data & Methods');
  await edgeCheck(page, tag + ' methods', '[data-view-scroll] h1');
  {
    const name = tag + ' methods count by the guidance\'s wording matches the release summary';
    if (GUIDANCE === null) skipped.push({ check: name, reason: GUIDANCE_MISSING });
    else {
      const def = await page.evaluate(() => { const dt = [...document.querySelectorAll('#m-scores dt')].find(d => d.textContent.trim() === 'VME and ME'); return dt?.nextElementSibling?.textContent.replace(/\s+/g, ' ') ?? null; });
      const m = (def ?? '').match(/Read by that wording, over the isolates each learner called, the best learner meets on ([\d,]+) of ([\d,]+) datasets? with lineages held out\./);
      const want = guidanceRow('best:lineage');
      check(name, m !== null && toCount(m[1]) === want.meets && toCount(m[2]) === want.n, { release: GUIDANCE.release, read: m === null ? null : [toCount(m[1]), toCount(m[2])], want: [want.meets, want.n] });
    }
  }
  await controlCheck(page, tag + ' methods');
  const tracks = await page.evaluate(() => [...document.querySelectorAll('#m-funnel td span.relative')].map(e => getComputedStyle(e).backgroundColor).filter(c => c !== 'rgba(0, 0, 0, 0)').length);
  check(tag + ' funnel bars sit on tracks', tracks >= 10, tracks);
  await scrollTo(page, '#m-funnel'); await shot(page, tag + '-definitions');
  check(tag + ' top bar stays in view after scrolling Methods', await page.evaluate(() => document.querySelector('header').getBoundingClientRect().top >= -1 && document.scrollingElement.scrollTop === 0), null);
  // Assays the report could not score: one entry each in Methods, or no such
  // part when there are none. An entry's lines are its name, its datasets
  // ("<species> vs <drug>", comma-separated, each matched exactly) and its
  // reason.
  const unscored = manifest.unsupported_panels ?? [];
  const unscoredPart = await page.evaluate(() => { const p = document.getElementById('m-unsupported'); return p ? [...p.querySelectorAll('li')].map(li => [...li.querySelectorAll('p')].map(x => x.textContent)) : null; });
  check(tag + ' assays not scored listed in Methods', unscored.length === 0 ? unscoredPart === null
    : unscoredPart !== null && unscoredPart.length === new Set(unscored.map(u => u.panel)).size && unscored.every(u => unscoredPart.some(([name = '', sets = '', ...rest]) => name === u.name && sets.split(', ').includes(datasetTitle(u)) && rest.join(' ').includes('Not scored.') && rest.join(' ').toLowerCase().includes(u.reason.toLowerCase()))), { unscored: unscored.map(u => [u.panel, u.dataset]), part: unscoredPart });
  if (unscoredPart !== null) { await scrollTo(page, '#m-unsupported'); await shot(page, tag + '-unscored'); }
  await page.goto(BASE);
  await page.locator('#corpus').waitFor();
  for (const button of await page.getByRole('button', { name: 'Method and coverage' }).all()) await button.click();
  await reachCheck(page, tag + ' overview with every Method and coverage open');
  // Each finding keeps its longer table one control away. With all of them
  // open, a wide table scrolls inside its own frame and the page keeps its
  // width.
  // A finding offers its table only when the tree has a dataset to count in
  // it, so the tree decides which controls the page owes.
  const longerTables = [
    ['Show best learner against each learner', manifest.datasets.some(d => d.best_lineage !== null)],
    ['Show each learner against the bars', expected.stats.barsCounted > 0],
    ['Show rules grouped by their own score', expected.catalog.length > 0],
  ];
  const tablesOpened = [];
  for (const [name] of longerTables) { const b = page.getByRole('button', { name, exact: true }); if (await b.count() === 1) { await b.click(); tablesOpened.push(name); } }
  const tablesOwed = longerTables.filter(([, owed]) => owed).map(([name]) => name);
  check(tag + ' overview longer tables one control away', tablesOpened.length === tablesOwed.length && tablesOwed.every(name => tablesOpened.includes(name)), { tablesOpened, tablesOwed });
  await reachCheck(page, tag + ' overview with every longer table open');
  const wideTables = await page.evaluate(() => [...document.querySelectorAll('[data-view-scroll] section[id] .scroll-frame [role="region"]')].map(r => ({ label: r.getAttribute('aria-label'), section: r.closest('section[id]').id, scroll: r.scrollWidth, client: r.clientWidth, table: r.querySelector('table')?.getBoundingClientRect().width ?? 0, columns: r.querySelectorAll('thead th').length })));
  // With all three tables open the page holds at least four tables of five
  // or more columns; a tree that owes fewer is held to the width alone.
  check(tag + ' overview tables of five or more columns keep a readable width', (tablesOwed.length < longerTables.length || wideTables.filter(t => t.columns > 4).length >= 4) && wideTables.filter(t => t.columns > 4).every(t => t.table >= 543), wideTables);
  for (const id of ['learners', 'bars', 'catalog']) { await scrollTo(page, '#' + id); await shot(page, tag + '-overview-' + id + '-open'); }
  // The bar counts are the page's own reading of the index: the three
  // outcomes add up to the datasets counted.
  const barTiles = await page.evaluate(() => [...document.querySelectorAll('#bars ol li')].map(li => ({ value: Number(li.firstElementChild.textContent.replace(/,/g, '')), unit: li.lastElementChild.textContent })));
  // A tree with no dataset to count draws no tiles.
  check(tag + ' overview bar outcomes add up', expected.stats.barsCounted === 0 ? barTiles.length === 0 : barTiles.length === 4 && barTiles[0].value === expected.stats.barsCounted && barTiles[1].value + barTiles[2].value + barTiles[3].value === barTiles[0].value, { barTiles, counted: expected.stats.barsCounted });
  await wordingCheck(page, tag);
  // The primer is a section of its own above the ten numbered ones, under
  // the same hairline.
  const rules = await page.evaluate(() => { const blocks = [...document.querySelectorAll('[data-view-scroll] section[id]')]; return blocks.map(b => ({ id: b.id, rule: getComputedStyle(b).borderTopStyle })); });
  const numbered = rules.filter(r => r.id !== 'primer');
  check(tag + ' one hairline between findings', numbered.length === 10 && rules.length === 11 && rules[0].id === 'primer' && rules.every(r => r.rule === 'solid'), rules);
  await page.goto(BASE);
  await page.locator('#corpus').waitFor();
  await headerCheck(page, tag + ' overview', 'What GRM found');
  await edgeCheck(page, tag + ' overview', '#corpus h2');
  const counts = await page.evaluate(() => [...document.querySelectorAll('#corpus ul li')].map(li => ({ text: li.firstElementChild.textContent, scroll: li.firstElementChild.scrollWidth, client: li.firstElementChild.clientWidth })));
  check(tag + ' overview corpus counts fit their cells', counts.length === 4 && counts.every(c => c.scroll <= c.client + 1), counts);
  const marks = await page.evaluate(() => { const c = document.querySelector('#lineage svg circle[data-dataset]'); return c ? getComputedStyle(c).fill : null; });
  const blue = await page.evaluate(() => { const d = document.createElement('div'); d.style.color = getComputedStyle(document.documentElement).getPropertyValue('--series-s'); document.body.append(d); const c = getComputedStyle(d).color; d.remove(); return c; });
  check(tag + ' chart marks in the blue series', marks === blue, { marks, blue });
  for (const id of ['lineage', 'learners', 'bars', 'catalog', 'panels', 'mic', 'outside', 'inputs']) { await scrollTo(page, '#' + id); await shot(page, tag + '-overview-' + id); }
  await outsideCheck(page, tag);
  await inputsCheck(page, tag);

  const micId = expected.mic[0].id;
  const micReport = JSON.parse(fs.readFileSync(path.join(FIXTURE, 'datasets', micId + '.json')));
  await page.goto(BASE + '/datasets?dataset=' + micId);
  const heat = page.getByRole('region', { name: 'Measured and predicted MIC' }); await heat.waitFor();
  await scrollTo(page, '#ds-mic'); await shot(page, tag + '-mic-heatmap');
  const cells = await heat.locator('svg g>rect').evaluateAll(es => es.map(e => ({ width: Number(e.getAttribute('width')), count: Number(e.nextElementSibling.textContent) })));
  check(tag + ' MIC cells and counts', cells.every(c => c.width >= 22) && JSON.stringify(cells.map(c => c.count)) === JSON.stringify(micReport.mic.confusion.map(c => c[2])), cells);
  const probeRow = manifest.panels.find(p => p.no_call !== undefined && (p.no_call.resistant + p.no_call.susceptible) > 0);
  if (probeRow !== undefined) {
    const probeReport = JSON.parse(fs.readFileSync(path.join(FIXTURE, 'datasets', probeRow.dataset + '.json')));
    const probe = probeReport.panels.find(p => p.id === probeRow.panel);
    await page.goto(BASE + '/datasets?dataset=' + probeRow.dataset);
    const noCall = page.getByText(/The probes give no call for/).first(); await noCall.waitFor(); await noCall.scrollIntoViewIfNeeded(); await page.waitForTimeout(150);
    check(tag + ' probe no-call details', await noCall.textContent().then(t => t.includes(`${probe.no_call.resistant} resistant and ${probe.no_call.susceptible} susceptible`)), await noCall.textContent());
    await shot(page, tag + '-probe-no-calls');
  } else skipped.push({ check: tag + ' probe no-call details', reason: 'no panel row with a probe no-call in this tree' });
}

async function states(browser, width, height, theme) {
  const tag = width + '-' + theme;
  const context = await browser.newContext({ viewport: { width, height }, reducedMotion: 'reduce' });
  await context.addInitScript(t => { if (/^https?:$/.test(location.protocol)) localStorage.setItem('grm-theme', t); }, theme);
  const page = await context.newPage();
  const cdp = await context.newCDPSession(page);
  await cdp.send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 1, mobile: false });
  page.on('pageerror', e => errors.push({ tag, state: true, message: e.message }));
  let release;
  const gate = new Promise(resolve => release = resolve);
  await page.route('**/data/index.json', async route => { await gate; await route.continue(); });
  await page.goto(BASE, { waitUntil: 'domcontentloaded' });
  await page.getByText('Loading the report', { exact: true }).waitFor();
  check(tag + ' loading aria-busy', await page.locator('[aria-busy="true"]').count() === 1, null);
  check(tag + ' loading runs the mark loop in the panel and the top bar', await page.locator('svg.grm-mark[data-motion="loop"]').count() === 2, await page.locator('svg.grm-mark').evaluateAll(es => es.map(e => e.getAttribute('data-motion'))));
  await shot(page, tag + '-loading'); release();
  await page.getByRole('heading', { name: 'What GRM found', exact: true }).waitFor();
  await page.unroute('**/data/index.json');
  await page.route('**/data/index.json', route => route.fulfill({ status: 500, body: 'unavailable' }));
  await page.reload(); await page.getByRole('alert').waitFor(); await shot(page, tag + '-manifest-error');
  await page.unroute('**/data/index.json');
  await page.route('**/data/datasets/' + inspectId + '.json', route => route.fulfill({ status: 404, body: 'not found' }));
  await page.goto(BASE + '/datasets?dataset=' + inspectId);
  await page.getByRole('alert').waitFor(); await shot(page, tag + '-dataset-error');
  await page.unroute('**/data/datasets/' + inspectId + '.json');
  const empty = { ...manifest, datasets: [], panels: [] };
  await page.route('**/data/index.json', route => route.fulfill({ json: empty }));
  await page.goto(BASE); await page.getByText('No dataset has a report yet', { exact: true }).waitFor(); await shot(page, tag + '-empty-manifest');
  await page.unroute('**/data/index.json');
  // A report from before the audit has no revision to make stale in one dataset.
  if (manifest.datasets[0].evaluation === undefined) {
    for (const name of [' mixed audit status is plain text without a hover wall', ' mixed audit row marker', ' mixed audit count and reason are visible']) skipped.push({ check: tag + name, reason: 'no dataset in this tree has an evaluation audit' });
  } else {
    const mixed = structuredClone(manifest); mixed.datasets[0].evaluation.revision = manifest.evaluation_revision - 1;
    await page.route('**/data/index.json', route => route.fulfill({ json: mixed }));
    await page.goto(BASE + '/datasets');
    const mixedStatus = page.getByText('Provisional metrics', { exact: true }).first();
    await mixedStatus.waitFor();
    await mixedStatus.hover();
    check(tag + ' mixed audit status is plain text without a hover wall', await mixedStatus.textContent() === 'Provisional metrics' && await page.getByRole('tooltip').count() === 0, await mixedStatus.textContent());
    check(tag + ' mixed audit row marker', await page.getByRole('region', { name: 'Datasets', exact: true }).getByText('Provisional', { exact: true }).count() === 1, null);
    await page.goto(BASE);
    const mixedNote = page.getByRole('note').first();
    await mixedNote.waitFor();
    const mixedNoteText = await mixedNote.textContent();
    check(tag + ' mixed audit count and reason are visible', mixedNoteText.includes('1 of ' + expected.datasets + ' dataset evaluations do not pass') && mixedNoteText.includes('current revision'), mixedNoteText);
    await shot(page, tag + '-mixed-audit');
    await page.unroute('**/data/index.json');
  }
  for (const [row, word] of stateRows) {
    await page.goto(BASE + '/datasets'); await page.getByRole('searchbox', { name: 'Search datasets' }).fill(row.antibiotic + ' ' + row.species);
    check(tag + ' ' + row.id + ' state ' + word, await page.getByText(word, { exact: true }).count() > 0, word);
    await shot(page, tag + '-state-' + row.id);
  }
  for (const word of ['Fold lacks a class', 'Too few lineages', 'No rule', 'No MIC']) if (!stateRows.some(([, w]) => w === word)) skipped.push({ check: tag + ' state ' + word, reason: 'no dataset in this tree is in this state' });
  await context.close();
}

async function present(page, tag, width, height) {
  await page.getByRole('button', { name: 'Present', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'GRM findings presentation' });
  await dialog.waitFor();
  check(tag + ' no deck chooser', !(await dialog.textContent()).includes('Choose a deck'), null);
  const slideSelect = dialog.getByRole('combobox', { name: /^Go to slide/ });
  // Webfonts and ResizeObserver measurements can add a continuation just
  // after Present opens. Read a stable outline before navigating it.
  await settle(page);
  await slideSelect.click();
  let slideCount = 0, previousCount = -1, stableCounts = 0;
  // Title, five primer steps, nine findings, the rule: sixteen slides, plus
  // the extra pages of any figure too dense for the viewport (the outline
  // labels each page "(k of n)"). Which figures continue depends on the
  // viewport and on the report's row counts, so the count is checked
  // against the labels the deck shows.
  const BASE_SLIDES = 16;
  for (let attempt = 0; attempt < 200; attempt++) {
    slideCount = await dialog.getByRole('option').count();
    if (slideCount === previousCount) stableCounts++;
    else { previousCount = slideCount; stableCounts = 0; }
    if (slideCount >= BASE_SLIDES && stableCounts >= 20) break;
    await page.waitForTimeout(50);
  }
  await page.keyboard.press('Escape');
  const labels = [];
  const panelSeen = [];
  const inputsSeen = [];
  const wordingSeen = [];
  for (let slide = 0; slide < slideCount; slide++) {
    await goSlide(page, dialog, slide);
    const title = await dialog.locator('main h2').first().textContent();
    if ((title ?? '').startsWith('Newly deposited genomes')) {
      const plots = await page.evaluate(outsidePlots);
      check(tag + ' slide ' + (slide + 1) + ' outside plots: no two axis labels touching', plots.length >= 2 && plots.every(p => p.touching.length === 0), plots);
    }
    { const w = ((await dialog.locator('main').textContent()) ?? '').replace(/\s+/g, ' ').match(/By the guidance's own wording, read after the results, it meets on ([\d,]+)\./); if (w !== null) wordingSeen.push(toCount(w[1])); }
    labels.push(await slideSelect.textContent());
    await shot(page, tag + '-slide-' + (slide + 1));
    check(tag + ' slide ' + (slide + 1) + ' selected in the outline', await slideNumber(dialog) === slide + 1, await slideNumber(dialog));
    const m = await page.evaluate(measure);
    check(tag + ' slide ' + (slide + 1) + ' true viewport', m.viewport[0] === width && m.viewport[1] === height, m.viewport);
    check(tag + ' slide ' + (slide + 1) + ' fits', m.main.scrollHeight <= m.main.height + 1 && m.main.scrollWidth <= m.main.width + 1 && m.section.scrollHeight <= m.section.height + 1, m);
    check(tag + ' slide ' + (slide + 1) + ' brand and title share an edge', m.titleX !== null && m.brandX !== null && Math.abs(m.titleX - m.brandX) <= 1, { titleX: m.titleX, brandX: m.brandX });
    // From 640px the footer is one row (Previous, outline, ticks, Next); on a
    // phone the outline and ticks take the first row and Previous and Next
    // the second. Either way no two controls overlap, all stay on screen, and
    // the outline shows its whole label.
    const foot = await page.evaluate(() => {
      const f = document.querySelector('dialog.presentation-dialog footer');
      const buttons = [...f.querySelectorAll('button')];
      const box = (name, e) => { if (!e) return null; const r = e.getBoundingClientRect(); const k = v => Math.round(v * 10) / 10; return { name, left: k(r.left), right: k(r.right), top: k(r.top), bottom: k(r.bottom), height: k(r.height) }; };
      const label = f.querySelector('[role="combobox"] [data-dropdown-label]');
      return {
        boxes: [box('previous', buttons.find(b => b.textContent === 'Previous')), box('outline', f.querySelector('[role="combobox"]')), box('ticks', f.querySelector('[data-beats]')), box('next', buttons.find(b => b.textContent === 'Next' || b.textContent === 'Finish'))].filter(Boolean),
        label: label ? { text: label.textContent, scroll: label.scrollWidth, client: label.clientWidth, overflow: getComputedStyle(label).textOverflow, white: getComputedStyle(label).whiteSpace } : null,
      };
    });
    const fb = foot.boxes, byName = n => fb.find(b => b.name === n);
    const apart = fb.every((a, i) => fb.every((b, j) => i === j || a.right <= b.left + 0.5 || b.right <= a.left + 0.5 || a.bottom <= b.top + 0.5 || b.bottom <= a.top + 0.5));
    const rows = width < 640 ? byName('outline').bottom <= Math.min(byName('previous').top, byName('next').top) + 0.5 && Math.abs(byName('previous').top - byName('next').top) <= 1 : fb.every((r, i) => i === 0 || fb[i - 1].right <= r.left + 0.5);
    check(tag + ' slide ' + (slide + 1) + ' footer controls keep their own space', apart && rows && byName('previous').right <= byName('next').left && fb.every(b => b.left >= -0.5 && b.right <= width + 0.5 && b.bottom <= height + 0.5), foot.boxes);
    check(tag + ' slide ' + (slide + 1) + ' outline label is whole', foot.label !== null && foot.label.scroll <= foot.label.client + 1 && foot.label.overflow !== 'ellipsis' && foot.label.white !== 'nowrap' && byName('outline').height >= 27.5, foot.label);
    const deck = await page.evaluate(() => { const t = document.querySelector('dialog.presentation-dialog header [data-deck-title]'); const h = t?.closest('header'); return t ? { text: t.textContent.trim(), scroll: t.scrollWidth, client: t.clientWidth, overflow: getComputedStyle(t).textOverflow, right: t.getBoundingClientRect().right, close: h.querySelector('button').getBoundingClientRect().left } : null; });
    check(tag + ' slide ' + (slide + 1) + ' deck title is whole', deck !== null && deck.text === 'What GRM found · Genomic Resistance Mapping' && deck.scroll <= deck.client + 1 && deck.overflow !== 'ellipsis' && deck.right <= deck.close, deck);
    // A fixed figure sits one 24px gap under its answer, not centred in the
    // room the slide leaves.
    const gap = await page.evaluate(() => { const s = document.querySelector('dialog.presentation-dialog main section'); const box = s?.querySelector('[data-fixed-box]'); if (!box) return null; const answer = s.querySelector('h2 + div p, h2 + p'); const fig = box.querySelector('[data-figure]'); return { answerBottom: answer.getBoundingClientRect().bottom, figureTop: fig.getBoundingClientRect().top, gap: fig.getBoundingClientRect().top - answer.getBoundingClientRect().bottom }; });
    if (gap !== null) check(tag + ' slide ' + (slide + 1) + ' fixed figure one gap under its answer', Math.abs(gap.gap - 24) <= 2, gap);
    if (slide === 3) await splitTracks(page, tag + ' slide 4', 'dialog.presentation-dialog main section');
    if (slide === 0) {
      const t = await page.evaluate(() => ({ mark: document.querySelector('dialog.presentation-dialog main svg.grm-mark')?.getAttribute('data-motion'), title: document.querySelector('dialog.presentation-dialog main h2')?.textContent, text: document.querySelector('dialog.presentation-dialog main')?.textContent }));
      check(tag + ' title slide: mark choreography, wordmark, date, no extra corpus count', t.mark === 'present' && t.title === 'Genomic Resistance Mapping' && t.text.includes('Report built') && !t.text.includes('BV-BRC lab records'), t);
      continue;
    }
    if (slide >= 1 && slide <= 5) {
      const primer = await page.evaluate(() => { const s = document.querySelector('dialog.presentation-dialog main section'); return { figure: s?.querySelector('[data-primer-figure]')?.getAttribute('data-primer-figure') ?? null, interactive: s ? s.querySelectorAll('a, button, input, [tabindex="0"]').length : -1 }; });
      const want = ['dataset', 'feature', 'splits', 'errors', 'scope'][slide - 1];
      check(tag + ' slide ' + (slide + 1) + ' primer step with its figure and no controls of its own', primer.figure === want && primer.interactive === 0, { ...primer, title });
    }
    check(tag + ' slide ' + (slide + 1) + ' figure above caption', m.contentBottom <= m.caption.y - 8 && (m.fit === null || m.fit.scrollHeight <= m.fit.height + 1), { contentBottom: m.contentBottom, caption: m.caption, fit: m.fit });
    check(tag + ' slide ' + (slide + 1) + ' labels inside', m.textOutside.length === 0, m.textOutside);
    const finding = title.includes('unseen lineages') ? 1 : title.includes('Known-gene rules') ? 2 : title.includes('PCR panel') ? 3 : title.includes('(MIC)') ? 4 : 0;
    // A paged panel figure draws its share of the rows on each page; the
    // pages together must draw every lead-panel dataset once.
    const paged = /\((\d+) of (\d+)\)/.exec(labels[slide] ?? '');
    const panelPer = paged === null ? expected.panels.length : Math.ceil(expected.panels.length / Number(paged[2]));
    const panelShown = paged === null ? expected.panels.length : Math.min(panelPer, expected.panels.length - (Number(paged[1]) - 1) * panelPer);
    const count = [0, expected.paired.length, expected.catalog.length, panelShown, expected.mic.length][finding] ?? 0;
    if (finding > 0) check(tag + ' slide ' + (slide + 1) + ' all dataset marks', m.marks === count, { actual: m.marks, expected: count, title, page: paged?.[0] ?? null });
    await referenceChecks(page, tag, finding, panelSeen, paged !== null);
    if (title === 'Genomes and lab calls') {
      const stats = await page.evaluate(() => [...document.querySelectorAll('dialog.presentation-dialog main ol li')].map(li => { const v = li.firstElementChild; const r = v.getBoundingClientRect(); return { text: v.textContent, left: r.left, right: r.right, top: r.top, bottom: r.bottom, cellRight: li.getBoundingClientRect().right, scroll: v.scrollWidth, client: v.clientWidth }; }));
      const apart = stats.every((a, i) => stats.every((b, j) => i === j || a.right <= b.left || b.right <= a.left || a.bottom <= b.top || b.bottom <= a.top));
      check(tag + ' corpus counts do not overlap', stats.length === 4 && apart && stats.every(s => s.right <= s.cellRight + 1 && s.scroll <= s.client + 1), stats);
    }
    if (title === 'Newly deposited genomes') {
      // Both populations on one slide, or one on each page of two, never
      // pooled: each plot draws one population's rows on one measure.
      const plots = await page.evaluate(() => [...document.querySelectorAll('dialog.presentation-dialog main section [role="listbox"]')].map(l => l.querySelectorAll('[data-outside]').length));
      const shown = paged === null ? outsideExpected : [outsideExpected[Number(paged[1]) - 1]];
      check(tag + ' slide ' + (slide + 1) + ' outside marks by population', JSON.stringify(plots) === JSON.stringify(shown.flatMap(p => [p.read, p.read])), { plots, page: paged?.[0] ?? null });
    }
    if (title === 'Input settings against the default') {
      // The rows in order, a share on each page; the pages together show
      // every row once, and no internal name.
      const shownRows = await page.evaluate(() => [...document.querySelectorAll('dialog.presentation-dialog main section tr[data-input]')].map(tr => tr.getAttribute('data-input')));
      const text = await page.evaluate(() => document.querySelector('dialog.presentation-dialog main section')?.textContent ?? '');
      inputsSeen.push(...shownRows);
      check(tag + ' slide ' + (slide + 1) + ' input rows and no internal name', shownRows.length > 0 && !INTERNAL.test(text), { shownRows, page: paged?.[0] ?? null, internal: text.match(INTERNAL)?.[0] ?? null });
    }
    const beats = await dialog.locator('[data-beats] span').evaluateAll(es => es.map(e => getComputedStyle(e).backgroundColor));
    const frameLabel = await slideSelect.textContent();
    const continuation = /\((\d+) of \d+\)/.exec(frameLabel ?? '');
    if (continuation !== null && Number(continuation[1]) > 1) {
      // A continuation stages only its figure; with one beat there are no
      // ticks to land, rather than a row of unfinished ticks.
      check(tag + ' slide ' + (slide + 1) + ' beat ticks all landed', beats.length === 0, beats);
    } else {
      check(tag + ' slide ' + (slide + 1) + ' beat ticks all landed', beats.length >= 2 && new Set(beats).size === 1, beats);
    }
    if (slide === 1) {
      const before = await slideNumber(dialog);
      await dialog.locator('[data-caption]').first().click();
      check(tag + ' caption click does not advance', await slideNumber(dialog) === before, { before, after: await slideNumber(dialog) });
    }
  }
  // Each paged figure is one run of adjacent slides labelled 1 to n.
  const pages = labels.map(l => /\((\d+) of (\d+)\)$/.exec(l ?? '')).filter(Boolean);
  const runs = [];
  labels.forEach((l, i) => {
    const p = /\((\d+) of (\d+)\)$/.exec(l ?? '');
    if (p === null) return;
    if (p[1] === '1' || runs.length === 0 || runs[runs.length - 1].end !== i - 1) runs.push({ end: i, pages: [] });
    runs[runs.length - 1].end = i;
    runs[runs.length - 1].pages.push([Number(p[1]), Number(p[2])]);
  });
  const wellFormed = runs.every(r => r.pages.length === r.pages[0][1] && r.pages.every(([k, n], i) => k === i + 1 && n === r.pages[0][1]));
  check(tag + ' expected slide count: sixteen slides plus each figure page past the first, labelled in order', slideCount === BASE_SLIDES + pages.filter(p => p[1] !== '1').length && wellFormed, { slideCount, labels });
  check(tag + ' panel pages together draw every lead-panel dataset once', panelSeen.length === expected.panels.length && JSON.stringify([...panelSeen].sort()) === JSON.stringify(expected.panels.map(p => p.dataset).sort()), { seen: panelSeen.length, expected: expected.panels.length });
  check(tag + ' input pages together show every row once, in order', JSON.stringify(inputsSeen) === JSON.stringify(inputsExpected.map(e => e.id)), { seen: inputsSeen });
  // The bars slide's caption gives the best learner's count by the
  // guidance's own wording with lineages held out, once.
  if (GUIDANCE === null) skipped.push({ check: tag + ' bars slide count by the guidance\'s wording matches the release summary', reason: GUIDANCE_MISSING });
  else check(tag + ' bars slide count by the guidance\'s wording matches the release summary', wordingSeen.length === 1 && wordingSeen[0] === guidanceRow('best:lineage').meets, { release: GUIDANCE.release, seen: wordingSeen, want: guidanceRow('best:lineage').meets });
  const closeButton = dialog.getByRole('button', { name: 'Close', exact: true });
  await closeButton.focus(); await page.keyboard.press('Home');
  check(tag + ' Home reaches the title slide', await slideNumber(dialog) === 1, await slideNumber(dialog));
  // Past the title and the five primer steps and the Data slide, to the
  // first chart.
  for (let i = 0; i < 7; i++) await page.keyboard.press('ArrowRight');
  const chart = dialog.getByRole('listbox'); await chart.focus(); await page.keyboard.press('ArrowRight');
  check(tag + ' chart owns arrow keys', await slideNumber(dialog) === 8 && await chart.getByRole('option').getAttribute('aria-posinset') === '2', await slideNumber(dialog));
  await shot(page, tag + '-chart-focus');
  await page.keyboard.press('Escape');
  check(tag + ' chart Escape keeps Present open', await dialog.count() === 1 && await dialog.getByRole('tooltip').count() === 0, null);
  await page.keyboard.press('Tab');
  check(tag + ' chart allows Tab exit', !await chart.evaluate(e => e.contains(document.activeElement)), await activeName(page));
  await closeButton.focus(); await page.keyboard.press('End');
  check(tag + ' End reaches the last slide', await slideNumber(dialog) === slideCount, await slideNumber(dialog));
  await page.keyboard.press('ArrowLeft');
  check(tag + ' ArrowLeft goes back', await slideNumber(dialog) === slideCount - 1, null);
  await page.keyboard.press('PageDown');
  check(tag + ' PageDown goes forward', await slideNumber(dialog) === slideCount, null);
  await page.keyboard.press('Escape');
  check(tag + ' Present focus return', (await activeName(page)) === 'Present', await activeName(page));
}

async function normalMotion(browser) {
  // The reduced-motion deck renders every slide's final text at once; read
  // slide 2's first, so the normal-motion deck can be waited on until its
  // staged beats land on that same text (up to STATE_MS), not for a fixed time.
  const reduced = await browser.newContext({ viewport: { width: 1440, height: 900 }, reducedMotion: 'reduce' });
  const page2 = await reduced.newPage();
  await page2.goto(BASE); await page2.getByRole('heading', { name: 'What GRM found', exact: true }).waitFor();
  await page2.getByRole('button', { name: 'Present', exact: true }).click();
  const d2 = page2.getByRole('dialog', { name: 'GRM findings presentation' }); await d2.waitFor();
  await d2.getByRole('button', { name: 'Close', exact: true }).focus(); await page2.keyboard.press('ArrowRight');
  const final = await d2.locator('main').innerText();
  const still = await page2.evaluate(() => document.getAnimations().filter(a => a.effect?.target?.closest?.('.grm-mark') && a.playState === 'running').length);
  check('reduced motion: the mark rests', still === 0, still);
  await reduced.close();
  const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, reducedMotion: 'no-preference' });
  const page = await context.newPage();
  await page.goto(BASE); await page.getByRole('heading', { name: 'What GRM found', exact: true }).waitFor();
  check('normal motion media preference', await page.evaluate(() => !matchMedia('(prefers-reduced-motion: reduce)').matches), null);
  // The title slide's beats land TITLE_REVEAL_MS after it opens (the title
  // delays in PresentMode.tsx). The first arrow must come before that to test
  // anything, so the check times the span from the deck's insertion (a
  // mutation observer) to the arrow's keydown; that span is never shorter
  // than the app's own. An arrow that moved on is a failure, unless the span
  // reached TITLE_REVEAL_MS: then the machine missed the reveal, and the
  // check closes the deck and opens it again, up to three times.
  const dialog = page.getByRole('dialog', { name: 'GRM findings presentation' });
  const arrows = [];
  for (;;) {
    await page.evaluate(() => {
      const t = window.__grmDeck = { open: null, key: null };
      const mo = new MutationObserver(() => { if (document.querySelector('dialog.presentation-dialog[open]')) { t.open = performance.now(); mo.disconnect(); } });
      mo.observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ['open'] });
      window.addEventListener('keydown', () => { if (t.key === null) t.key = performance.now(); }, { capture: true, once: true });
    });
    await page.getByRole('button', { name: 'Present', exact: true }).click();
    await dialog.waitFor();
    await dialog.getByRole('button', { name: 'Close', exact: true }).focus();
    await page.keyboard.press('ArrowRight');
    const t = await page.evaluate(() => window.__grmDeck);
    const a = { slide: await slideNumber(dialog), ms: t.open !== null && t.key !== null ? Math.round(t.key - t.open) : null };
    arrows.push(a);
    if (a.slide === 1 || a.ms === null || a.ms < TITLE_REVEAL_MS || arrows.length === 3) break;
    await page.keyboard.press('Escape');
    await dialog.waitFor({ state: 'detached' });
  }
  check('normal motion: an arrow mid-reveal lands the beats first', arrows.at(-1).slide === 1, { revealMs: TITLE_REVEAL_MS, arrows });
  await page.keyboard.press('ArrowRight');
  check('normal motion: the next arrow moves on', await slideNumber(dialog) === 2, await slideNumber(dialog));
  const t0 = Date.now();
  await page.waitForFunction(want => document.querySelector('dialog.presentation-dialog main')?.innerText === want, final, { timeout: STATE_MS }).catch(() => {});
  const landed = await dialog.locator('main').innerText();
  check('normal motion: staged beats land on the reduced-motion final text', landed === final, { ms: Date.now() - t0 });
  await page.keyboard.press('Escape');
  await context.close();
}

/** The two results slides on a short wide screen (1280×560): on every page
 *  of "Newly deposited genomes" and "Input settings against the default", no
 *  text box overlaps another, no figure is drawn over text, and nothing runs
 *  past the slide. The older slides that also crowd at this size are not
 *  checked here. */
/** On a slide, the outside study's live plots (a measuring copy has no
 *  keyboard layer): each plot's size, and any two of its axis labels that
 *  touch. */
const outsidePlots = () => [...document.querySelectorAll('dialog.presentation-dialog main [data-figure-fit] [data-outside-plot] svg.block')].filter(s => s.closest('[role="listbox"]') !== null).map(svg => {
  const t = [...svg.querySelectorAll('text')].filter(x => x.textContent !== '').map(x => ({ s: x.textContent, b: x.getBoundingClientRect() }));
  const touching = [];
  for (let i = 0; i < t.length; i++) for (let j = i + 1; j < t.length; j++) { const a = t[i].b, c = t[j].b; if (a.left < c.right && c.left < a.right && a.top < c.bottom && c.top < a.bottom) touching.push([t[i].s, t[j].s]); }
  return { side: Math.round(svg.getBoundingClientRect().width), labels: t.map(x => x.s), touching };
});
// The plot size the outside slide reaches at 1280x560 (114 px when this was
// set, with the hollow-dot sentence off the slide's caption and the
// caption's two qualifiers moved to the legend), less a small margin for
// font metrics.
const OUTSIDE_SHORT_MIN = 106;

async function shortSlides(browser) {
  const SHORT = [1280, 560];
  const TITLES = ['Newly deposited genomes', 'Input settings against the default'];
  const context = await browser.newContext({ viewport: { width: SHORT[0], height: SHORT[1] }, colorScheme: 'light', reducedMotion: 'reduce' });
  const page = await context.newPage();
  page.on('pageerror', error => errors.push({ width: SHORT[0], height: SHORT[1], message: error.message }));
  try {
    await page.goto(BASE); await page.getByRole('heading', { name: 'What GRM found', exact: true }).waitFor();
    await page.getByRole('button', { name: 'Present', exact: true }).click();
    const dialog = page.getByRole('dialog', { name: 'GRM findings presentation' }); await dialog.waitFor();
    await settle(page);
    // The outline grows as figures measure their pages; read it once stable,
    // with the slide chooser open (its options exist only then).
    await dialog.getByRole('combobox', { name: /^Go to slide/ }).click();
    let count = -1;
    for (let attempt = 0, stable = 0; attempt < 200 && stable < 20; attempt++) { const n = await dialog.getByRole('option').count(); stable = n === count ? stable + 1 : 0; count = n; await page.waitForTimeout(50); }
    await page.keyboard.press('Escape');
    const seen = Object.fromEntries(TITLES.map(t => [t, 0]));
    for (let slide = 0; slide < count; slide++) {
      await goSlide(page, dialog, slide);
      const r = await page.evaluate(titles => {
        const main = document.querySelector('dialog.presentation-dialog main');
        const title = titles.find(t => (main?.querySelector('h2')?.textContent ?? '').startsWith(t));
        if (!title) return { skip: true };
        const leaf = [...main.querySelectorAll('p, h1, h2, h3, span, th, td, li, figcaption')].filter(e => e.textContent.trim() !== '' && [...e.children].every(c => c.tagName.toLowerCase() === 'svg' || c.tagName === 'BR') && !e.closest('[inert], [aria-hidden="true"]') && getComputedStyle(e).visibility !== 'hidden');
        const boxes = leaf.map(e => ({ e, t: e.textContent.trim().slice(0, 40), b: e.getBoundingClientRect() })).filter(x => x.b.width > 0 && x.b.height > 0);
        const hit = (a, c) => a.left < c.right - 1 && c.left < a.right - 1 && a.top < c.bottom - 1 && c.top < a.bottom - 1;
        const over = [];
        for (let i = 0; i < boxes.length; i++) for (let j = i + 1; j < boxes.length; j++) if (hit(boxes[i].b, boxes[j].b) && !boxes[i].e.contains(boxes[j].e) && !boxes[j].e.contains(boxes[i].e)) over.push([boxes[i].t, boxes[j].t]);
        for (const g of main.querySelectorAll('svg')) {
          const a = g.getBoundingClientRect();
          if (g.closest('[inert]') || a.width <= 20 || a.height <= 20) continue;
          for (const x of boxes) if (!x.e.contains(g) && hit(a, x.b)) over.push(['<svg>', x.t]);
        }
        const mb = main.getBoundingClientRect();
        return { title, over, past: boxes.filter(x => x.b.bottom > mb.bottom + 1 || x.b.right > mb.right + 1).map(x => x.t) };
      }, TITLES);
      if (r.skip) continue;
      seen[r.title]++;
      check('1280x560 slide ' + (slide + 1) + ' (' + r.title + '): no text over text or figure, nothing past the slide', r.over.length === 0 && r.past.length === 0, { over: r.over.slice(0, 4), past: r.past.slice(0, 4) });
      if (r.title === TITLES[0]) {
        const plots = await page.evaluate(outsidePlots);
        check('1280x560 slide ' + (slide + 1) + ' (' + r.title + '): every plot at least ' + OUTSIDE_SHORT_MIN + ' px, no two axis labels touching', plots.length >= 2 && plots.every(p => p.side >= OUTSIDE_SHORT_MIN && p.touching.length === 0), { min: OUTSIDE_SHORT_MIN, plots });
      }
    }
    check('1280x560 both results slides reached', TITLES.every(t => seen[t] > 0), seen);
  } catch (error) { errors.push({ width: SHORT[0], height: SHORT[1], message: error.stack }); }
  await context.close();
}

async function boundaries(browser) {
  const context = await browser.newContext({ viewport: { width: 900, height: 900 }, reducedMotion: 'reduce', hasTouch: true });
  await context.addInitScript(() => { if (/^https?:$/.test(location.protocol)) localStorage.setItem('grm-filter-rail-open', '1'); });
  const page = await context.newPage();
  for (const width of [899, 900, 1279, 1280, 1439, 1440, 1680]) {
    await page.setViewportSize({ width, height: 900 }); await page.goto(BASE + '/datasets'); await page.getByRole('searchbox', { name: 'Search datasets' }).waitFor();
    const docked = await page.locator('aside[data-view-options-mode="inline"]').count();
    check('rail at ' + width, docked === (width >= 1440 ? 1 : 0), docked);
    const bar = await page.evaluate(() => { const nav = document.querySelector('nav[aria-label="Views"], nav[aria-label="View"]'); const brand = document.querySelector('header a[href="/"]'); return { pills: document.querySelector('nav[aria-label="Views"]') !== null, navMid: nav.getBoundingClientRect().top + nav.getBoundingClientRect().height / 2, brandMid: brand.getBoundingClientRect().top + brand.getBoundingClientRect().height / 2, barHeight: document.querySelector('header').getBoundingClientRect().height }; });
    check('top bar one row at ' + width, Math.abs(bar.navMid - bar.brandMid) <= 1.5 && bar.barHeight === 56, bar);
    await shot(page, 'boundary-' + width);
  }
  for (const [width, zoom] of [[720, 200], [360, 400]]) {
    await page.setViewportSize({ width, height: 900 }); await page.goto(BASE + '/datasets'); await page.getByRole('searchbox', { name: 'Search datasets' }).waitFor();
    check(zoom + '% equivalent reflow', await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), await page.evaluate(() => ({ width: innerWidth, scroll: document.documentElement.scrollWidth })));
    await shot(page, 'reflow-' + zoom);
  }
  await page.setViewportSize({ width: 390, height: 844 }); await page.emulateMedia({ forcedColors: 'active' }); await page.goto(BASE + '/datasets'); await page.getByRole('searchbox', { name: 'Search datasets' }).waitFor();
  const targets = await page.evaluate(() => [...document.querySelectorAll('button,input[type="search"]')].filter(e => !e.disabled && e.getBoundingClientRect().width > 0 && e.getBoundingClientRect().top < innerHeight && e.getBoundingClientRect().bottom > 0).map(e => ({ name: e.getAttribute('aria-label') ?? e.textContent, width: e.getBoundingClientRect().width, height: e.getBoundingClientRect().height })).filter(e => e.width < 44 || e.height < 44));
  check('coarse pointer 44px targets', targets.length === 0, targets);
  await page.getByRole('button', { name: 'Random', exact: true }).focus(); await shot(page, 'forced-colors-focus');
  check('forced colors: judged values keep plain readable text and no ring', await page.locator('[data-verdict]').first().evaluate(e => getComputedStyle(e).boxShadow === 'none' && e.textContent.trim().length > 0), null);
  const forcedMarks = await page.evaluate(() => [...document.querySelectorAll('.dataset-sub svg.class-glyph')].map(g => { const s = g.querySelector('circle, rect, path'); return { shape: g.getAttribute('data-shape'), paint: g.getAttribute('data-shape') === 'ring' ? getComputedStyle(s).stroke : getComputedStyle(s).fill }; }));
  check('forced colors: class marks share one system ink and keep distinct shapes', forcedMarks.length > 0 && new Set(forcedMarks.map(m => m.paint)).size === 1 && new Set(forcedMarks.map(m => m.shape)).size >= Math.min(2, new Set(manifest.datasets.map(d => d.drug_class)).size), forcedMarks.slice(0, 6));
  await longRuleCheck(page, 'touch 390', ['tap']);
  await holdCheck(page, context, 'touch 390');
  await page.goto(BASE);
  await page.getByRole('button', { name: 'Present', exact: true }).click();
  const deck = page.getByRole('dialog', { name: 'GRM findings presentation' });
  await deck.waitFor();
  await deck.getByRole('combobox', { name: /^Go to slide/ }).waitFor();
  await settle(page);
  for (const [slide, selector, name] of [[3, '[data-hidden], [data-hidden-swatch]', 'splits'], [4, '[data-error]', 'errors']]) {
    await goSlide(page, deck, slide);
    const marks = await page.evaluate((selector) => [...document.querySelectorAll('dialog.presentation-dialog main ' + selector)].map(e => ({ adjust: getComputedStyle(e).forcedColorAdjust, fill: getComputedStyle(e).backgroundColor })), selector);
    check('forced colors: ' + name + ' semantic fills stay visible', marks.length > 0 && marks.every(m => m.adjust === 'none' && m.fill !== 'rgba(0, 0, 0, 0)' && m.fill !== 'transparent'), marks);
    await shot(page, 'forced-colors-primer-' + name);
  }
  await context.close();
}

(async () => {
  const browser = await chromium.launch({ executablePath: process.env.GRM_CHROMIUM || undefined, headless: true, args: ['--no-sandbox', '--disable-dev-shm-usage'] });
  try {
    const favicon = await fetch(BASE + '/favicon.svg');
    check('favicon is the mark', favicon.ok && (favicon.headers.get('content-type') ?? '').includes('svg') && (await favicon.text()).includes('#ff8c66'), favicon.status);
    for (const [width, height] of [[1280, 720], [1440, 900], [390, 844]]) for (const theme of ['light', 'dark']) {
      const context = await browser.newContext({ viewport: { width, height }, colorScheme: theme, reducedMotion: 'reduce' });
      // Init scripts also run on each new page's first about:blank, which has
      // no storage, so touching localStorage there is a page error. Every
      // init script here sets it only on the app's own http(s) pages.
      await context.addInitScript(theme => { if (!/^https?:$/.test(location.protocol)) return; localStorage.setItem('grm-theme', theme); localStorage.setItem('grm-filter-rail-open', '1'); }, theme);
      const page = await context.newPage();
      const cdp = await context.newCDPSession(page);
      await cdp.send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 1, mobile: false });
      page.on('pageerror', error => errors.push({ theme, width, message: error.message }));
      page.on('response', r => { if (r.status() >= 400 && !r.url().endsWith('/data/research.json') && !(r.status() === 404 && r.url().includes('/data/intervals/'))) errors.push({ theme, width, url: r.url(), status: r.status() }); });
      const tag = width + '-' + theme;
      try {
        await page.goto(BASE);
        await page.getByRole('heading', { name: 'What GRM found', exact: true }).waitFor();
        await present(page, tag, width, height);
        if (width === 1280) { await interactions(page, tag, width); await widestIntervalFit(page, tag); await stripFit(page, tag); await panelFit(page, tag, width); await context.close(); continue; }
        for (const [route, title, key] of [['/', 'What GRM found', 'overview'], ['/datasets', 'Datasets', 'datasets'], ['/species/' + speciesSlug, speciesName, 'species'], ['/methods', 'Data & Methods', 'methods']]) {
          await page.goto(BASE + route);
          await page.getByRole('heading', { name: title, exact: true }).first().waitFor();
          if (route.startsWith('/species')) await page.getByRole('grid', { name: /lineage block/ }).first().waitFor();
          if (route === '/methods') await page.getByRole('heading', { name: /Models and evaluation/ }).waitFor();
          await shot(page, tag + '-' + key);
          if (key === 'overview') await primerCheck(page, tag);
          await reachCheck(page, tag + ' ' + key);
          check(tag + ' ' + key + ' document has no horizontal overflow', await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), await page.evaluate(() => ({ width: innerWidth, scroll: document.documentElement.scrollWidth })));
          await tallShot(page, tag + '-' + key + '-full', width, height);
          check(tag + ' ' + key + ' page overflow', await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), await page.evaluate(() => ({ width: innerWidth, scroll: document.documentElement.scrollWidth })));
          await headerCheck(page, tag + ' ' + key, title);
          check(tag + ' ' + key + ' 1560px content rail', await page.evaluate(() => { const shell = document.querySelector('[data-view-scroll] > div'); return getComputedStyle(shell).maxWidth === '1560px'; }), null);
        }
        if (theme === 'light') await smallMaps(page, tag);
        await page.goto(BASE + '/datasets?dataset=' + inspectId);
        await page.locator('#ds-learners').waitFor();
        await shot(page, tag + '-inspector');
        await intervalCheck(page, tag);
        await interactions(page, tag, width);
        if (width >= 1280) { await widestIntervalFit(page, tag); await stripFit(page, tag); }
        await panelFit(page, tag, width);
        await states(browser, width, height, theme);
      } catch (error) { errors.push({ theme, width, message: error.stack }); }
      await context.close();
    }
    await boundaries(browser);
    await shortSlides(browser);
    await normalMotion(browser);
  } finally { await browser.close(); }
  const receipt = { source: { fixture: FIXTURE, schema: manifest.schema, generated: manifest.generated_at, revision: manifest.evaluation_revision, datasets: expected.datasets, panelRows: manifest.panels.length }, independent: expected.stats, stateRows: stateRows.map(([d, w]) => [d.id, w]), skipped, checks, screenshots, errors, passed: checks.every(c => c.pass) && errors.length === 0 };
  fs.writeFileSync(path.join(OUT, 'receipt.json'), JSON.stringify(receipt, null, 2));
  console.log(JSON.stringify({ passed: receipt.passed, checks: checks.length, skipped, failures: checks.filter(c => !c.pass).map(c => ({ name: c.name, actual: c.actual })), errors, screenshots: screenshots.length }, null, 2));
  process.exitCode = receipt.passed ? 0 : 1;
})().catch(e => { console.error(e); process.exitCode = 1; });
