// Run on Vast against an already built artifact. The report is read-only.
// node browser-check.cjs <playwright-core> <output-dir> <base-url> <fixture-dir>
const fs = require('node:fs');
const path = require('node:path');
const { chromium } = require(process.argv[2]);
const OUT = process.argv[3], BASE = process.argv[4], FIXTURE = process.argv[5];
fs.mkdirSync(OUT, { recursive: true });
const manifest = JSON.parse(fs.readFileSync(path.join(FIXTURE, 'index.json')));
const checks = [], screenshots = [], errors = [];
const check = (name, pass, actual) => checks.push({ name, pass: Boolean(pass), actual });
const median = values => { const s=[...values].sort((a,b)=>a-b), k=s.length>>1;return s.length%2?s[k]:(s[k-1]+s[k])/2; };
const expected = {
  datasets: manifest.datasets.length,
  paired: manifest.datasets.filter(d => d.best_random !== null && d.best_lineage !== null),
  catalog: manifest.datasets.filter(d => d.catalog !== null && d.best_lineage !== null),
  mic: manifest.datasets.filter(d => d.mic_essential_agreement !== null),
  panels: manifest.panels.filter(p => p.panel === 'carbapenemase'),
};
const gaps=expected.paired.map(d=>d.best_random-d.best_lineage);
expected.stats={paired:expected.paired.length,medianRandom:median(expected.paired.map(d=>d.best_random)),medianLineage:median(expected.paired.map(d=>d.best_lineage)),medianGap:median(gaps),above5:gaps.filter(v=>v>0.05).length,catalogAhead:expected.catalog.filter(d=>d.catalog.balanced_accuracy-d.best_lineage>0.02).length,catalogWithin:expected.catalog.filter(d=>Math.abs(d.catalog.balanced_accuracy-d.best_lineage)<=0.02).length,catalogBehind:expected.catalog.filter(d=>d.catalog.balanced_accuracy-d.best_lineage< -0.02).length,micAbove90:expected.mic.filter(d=>d.mic_essential_agreement>=0.9).length};
const settle = async page => {
  await page.evaluate(() => document.fonts.ready);
  await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
};
const shot = async (page, name) => {
  await settle(page);
  const file = path.join(OUT, name + '.png');
  await page.screenshot({ path: file });
  screenshots.push(file);
  console.log('Screenshot '+name);
};
const measure = () => {
  const main = document.querySelector('dialog.presentation-dialog main');
  const section = main?.querySelector('section');
  const rect = e => e ? ({ x:e.getBoundingClientRect().x, y:e.getBoundingClientRect().y, width:e.clientWidth, height:e.clientHeight, scrollWidth:e.scrollWidth, scrollHeight:e.scrollHeight, bottom:e.getBoundingClientRect().bottom }) : null;
  const texts = [...(section?.querySelectorAll('svg text') ?? [])].map(e=>({text:e.textContent, box:e.getBoundingClientRect().toJSON(), svg:e.closest('svg').getBoundingClientRect().toJSON()}));
  const textOutside = texts.filter(t => t.box.left < t.svg.left - 1 || t.box.right > t.svg.right + 1 || t.box.top < t.svg.top - 1 || t.box.bottom > t.svg.bottom + 1);
  const fit = section?.querySelector('[data-figure-fit]');
  const contentBottom = Math.max(0,...[...(fit?.querySelectorAll('*') ?? [])].filter(e=>!e.closest('[role="tooltip"]')&&!e.classList.contains('sr-only')).map(e=>e.getBoundingClientRect().bottom));
  return { viewport:[innerWidth,innerHeight], page:[document.documentElement.scrollWidth,document.documentElement.clientWidth], main:rect(main), section:rect(section), fit:rect(fit), contentBottom, caption:rect(section?.lastElementChild), textOutside, marks:section?.querySelectorAll('[data-dataset]').length ?? 0 };
};

async function referenceChecks(page, tag, slide) {
  if(![1,2,3,4].includes(slide))return;
  const geometry=await page.locator('dialog.presentation-dialog main svg').first().evaluate(svg=>({width:svg.getBoundingClientRect().width,height:svg.getBoundingClientRect().height,refs:[...svg.querySelectorAll('line[aria-hidden="true"]')].map(l=>({x:Number(l.getAttribute('x1')),y:Number(l.getAttribute('y1'))})),ids:[...svg.querySelectorAll('[data-dataset]')].map(e=>e.getAttribute('data-dataset')),text:svg.textContent}));
  const close=(actual,want)=>Math.abs(actual-want)<=1.5;
  let wanted=[];
  if(slide===1){
    const lo=Math.floor(Math.min(-0.02,...gaps)*100)/100,hi=Math.ceil(Math.max(0.1,...gaps)*100)/100;
    wanted=[expected.stats.medianGap,0.05,0].map(v=>20+(v-lo)/(hi-lo)*(geometry.width-40));
    check(tag+' lineage region labels',geometry.text.includes((gaps.length-expected.stats.above5)+' at or below 5')&&geometry.text.includes(expected.stats.above5+' above 5'),geometry.text);
  }
  if(slide===4){
    const values=expected.mic.map(d=>d.mic_essential_agreement),lo=Math.min(0.8,Math.floor(Math.min(...values)*10)/10);
    wanted=[0.9,median(values)].map(v=>20+(v-lo)/(1-lo)*(geometry.width-40));
  }
  if(slide===3){
    const counts=expected.panels.map(p=>({panel:p.counts,model:p.best_model_counts}));
    const panel=counts.reduce((a,c)=>a+c.panel.tp,0)/counts.reduce((a,c)=>a+c.panel.tp+c.panel.fn,0);
    const model=counts.reduce((a,c)=>a+c.model.tp,0)/counts.reduce((a,c)=>a+c.model.tp+c.model.fn,0);
    const start=Math.min(220,Math.round(geometry.width*0.5))+10;
    wanted=[panel,model].map(v=>start+v*(geometry.width-20-start));
    check(tag+' lead panel datasets',JSON.stringify([...geometry.ids].sort())===JSON.stringify(expected.panels.map(p=>p.dataset).sort()),geometry.ids);
  }
  if(slide===2){
    const lo=Math.max(0,Math.floor(Math.min(0.5,...expected.catalog.flatMap(d=>[d.catalog.balanced_accuracy,d.best_lineage]))*10)/10);
    check(tag+' catalog 50% reference',geometry.refs.some(r=>close(r.y,22+(1-0.5)/(1-lo)*(geometry.height-56))),geometry.refs);
  } else check(tag+' slide '+(slide+1)+' reference positions',wanted.length===geometry.refs.length&&wanted.every((x,i)=>close(geometry.refs[i].x,x)),{wanted,actual:geometry.refs});
}

async function contrastCheck(page, tag) {
  const values=await page.locator('[data-verdict]').evaluateAll(cells=>{
    const rgb=s=>(s.match(/[\d.]+/g)??[]).slice(0,3).map(Number);
    const luminance=rgb=>rgb.map(x=>{const v=x/255;return v<=0.04045?v/12.92:((v+0.055)/1.055)**2.4;}).reduce((s,v,i)=>s+v*[0.2126,0.7152,0.0722][i],0);
    return [...new Map(cells.map(e=>{const s=getComputedStyle(e),a=luminance(rgb(s.color)),b=luminance(rgb(s.backgroundColor));return [e.getAttribute('data-verdict'),{color:s.color,background:s.backgroundColor,ratio:(Math.max(a,b)+0.05)/(Math.min(a,b)+0.05)}];})).entries()];
  });
  check(tag+' status text contrast',values.length===2&&values.every(([,v])=>v.ratio>=4.5),values);
}

const scrollTo = async (page, selector) => {
  await page.locator(selector).evaluate(e => e.scrollIntoView({block:'start'}));
  await settle(page);
};
const activeName = page => page.evaluate(() => document.activeElement?.getAttribute('aria-label') ?? document.activeElement?.textContent?.trim());

async function interactions(page, tag, width) {
  await page.goto(BASE+'/datasets');
  const search=page.getByRole('searchbox',{name:'Search datasets'});
  await search.waitFor();
  await page.getByRole('button',{name:'GRM overview'}).focus();
  const focus=[];
  for(let i=0;i<9;i++){focus.push(await activeName(page));await page.keyboard.press('Tab');}
  check(tag+' top bar keyboard order',JSON.stringify(focus.slice(0,8))===JSON.stringify(['GRM overview','Overview','Datasets','Species','Methods','Report','Present',tag.includes('dark')?'Switch to light theme':'Switch to dark theme']),focus);
  await page.getByRole('button',{name:'Report',exact:true}).focus();
  await page.keyboard.press('Enter');
  await page.getByRole('dialog',{name:'Report details'}).waitFor();
  await shot(page,tag+'-report-dialog');
  await page.keyboard.press('Escape');
  check(tag+' report focus return',(await activeName(page))==='Report',await activeName(page));
  if(width<1440)await page.getByRole('button',{name:'Filters',exact:true}).click();
  await page.getByRole('checkbox',{name:/Has MIC regression/}).check();
  await shot(page,tag+'-filters');
  if(width<1440){await page.keyboard.press('Escape');check(tag+' filter focus return',(await activeName(page)).startsWith('Filters'),await activeName(page));}
  check(tag+' MIC filter count',await page.getByText(expected.mic.length+' of '+expected.datasets+' datasets',{exact:true}).count()===1,null);
  if(width<1440)await page.getByRole('button',{name:/^Filters/}).click();
  await page.getByRole('button',{name:'Clear',exact:true}).click();
  if(width<1440)await page.keyboard.press('Escape');
  await search.fill('zz-no-dataset-matches-zz');
  await page.getByRole('button',{name:'Clear filters and search'}).waitFor();
  await shot(page,tag+'-empty-filter');
  await page.getByRole('button',{name:'Clear filters and search'}).click();
  check(tag+' clear filter recovery',await search.inputValue()==='',null);
  await page.getByRole('button',{name:'Show all',exact:true}).click();
  const status=await page.getByRole('region',{name:'Datasets',exact:true}).evaluate(region=>{
    const rows=[...region.querySelectorAll('tbody tr')];
    return {rows:rows.length,vme:rows.filter(r=>r.children[5]?.getAttribute('data-verdict')==='meets').length,me:rows.filter(r=>r.children[6]?.getAttribute('data-verdict')==='meets').length,mic:rows.filter(r=>r.children[8]?.getAttribute('data-verdict')==='meets').length,zero:[...region.querySelectorAll('[data-verdict]')].some(e=>/0\.0%/.test(e.textContent)&&e.textContent.includes('meets')),headers:[...region.querySelectorAll('thead th')].map(e=>e.textContent),scope:[...region.querySelectorAll('thead th')].every(e=>e.scope==='col')};
  });
  const best=d=>d.learners.find(l=>l.learner===d.best_learner)?.random;
  check(tag+' status thresholds and counts',status.rows===expected.datasets&&status.vme===manifest.datasets.filter(d=>best(d)?.vme<=0.015).length&&status.me===manifest.datasets.filter(d=>best(d)?.me<=0.03).length&&status.mic===expected.mic.filter(d=>d.mic_essential_agreement>=0.9).length&&status.zero&&status.scope,status);
  await contrastCheck(page,tag);
  const table=page.getByRole('region',{name:'Datasets',exact:true});
  const left=await table.locator('tbody tr').first().locator('td').first().evaluate(e=>e.getBoundingClientRect().left);
  await table.evaluate(e=>e.scrollLeft=500);
  const shifted=await table.locator('tbody tr').first().locator('td').first().evaluate(e=>e.getBoundingClientRect().left);
  check(tag+' sticky dataset column',Math.abs(left-shifted)<=1.5,{left,shifted});
  await shot(page,tag+'-dataset-status');
  await table.evaluate(e=>e.scrollLeft=0);
  await page.getByRole('button',{name:'Random',exact:true}).focus();
  await page.keyboard.press('Enter');
  check(tag+' keyboard sorting',await table.locator('th[aria-sort="descending"]').textContent()==='Random',null);
  const rows=table.locator('[data-row-open]');
  await rows.nth(0).focus();await page.keyboard.press('ArrowDown');
  check(tag+' keyboard row movement',await rows.nth(1).evaluate(e=>e===document.activeElement),null);
  await shot(page,tag+'-row-focus');
  await page.keyboard.press('Enter');
  await page.locator('#ds-learners').waitFor();
  await page.keyboard.press('Escape');
  await page.waitForFunction(()=>!document.querySelector('dialog[open]'));
  check(tag+' inspector focus return',await rows.nth(1).evaluate(e=>e===document.activeElement),await activeName(page));

  await page.goto(BASE+'/datasets?dataset=escherichia_coli__ceftriaxone');
  await page.locator('#ds-model').waitFor();
  for(const id of ['ds-model','ds-mechanisms','ds-mic','ds-lineages']){
    await scrollTo(page,'#'+id);await shot(page,tag+'-'+id);
  }
  await page.getByRole('button',{name:'Method and coverage',exact:true}).click();
  await page.getByText('Evaluation audit',{exact:true}).last().scrollIntoViewIfNeeded();
  await shot(page,tag+'-method-fold');
  const owners=await page.getByRole('dialog').evaluate(d=>[...d.querySelectorAll('*')].filter(e=>getComputedStyle(e).overflowY==='auto'&&e.scrollHeight>e.clientHeight+1).map(e=>e.getAttribute('data-modal-scroll-owner')));
  check(tag+' inspector vertical scroll owner',owners.length===1&&owners[0]==='true',owners);
  await page.keyboard.press('Escape');
  await page.goto(BASE+'/species/escherichia_coli');
  await page.getByRole('grid').waitFor();
  await scrollTo(page,'#sp-cross');
  await shot(page,tag+'-cross-matrix');
  await page.getByRole('button',{name:'Table',exact:true}).click();
  await shot(page,tag+'-cross-table');
  check(tag+' cross-resistance table',await page.getByRole('region',{name:'Cross-resistance pairs'}).count()===1,null);
  await page.goto(BASE+'/tests');
  await page.getByRole('region',{name:'Panel results'}).waitFor();
  check(tag+' legacy tests redirect',new URL(page.url()).pathname==='/'&&new URL(page.url()).hash==='#panels',page.url());
  await shot(page,tag+'-panel-evidence');
  await page.goto(BASE+'/methods');
  await page.getByRole('heading',{name:'Antibiotic name resolution'}).waitFor();
  await scrollTo(page,'#m-design');await shot(page,tag+'-definitions');
  await scrollTo(page,'#m-unsupported');await shot(page,tag+'-unscored');
  // Every overview figure is also inspected in its page layout.
  await page.goto(BASE);
  for(const id of ['lineage','catalog','panels','mic']){await scrollTo(page,'#'+id);await shot(page,tag+'-overview-'+id);}

  // Exercise a populated MIC heatmap as well as the unavailable-model state.
  const micId=expected.mic[0].id;
  const micReport=JSON.parse(fs.readFileSync(path.join(FIXTURE,'datasets',micId+'.json')));
  await page.goto(BASE+'/datasets?dataset='+micId);
  const heat=page.getByRole('region',{name:'Measured and predicted MIC'});await heat.waitFor();
  await scrollTo(page,'#ds-mic');await shot(page,tag+'-mic-heatmap');
  const cells=await heat.locator('svg g>rect').evaluateAll(es=>es.map(e=>({width:Number(e.getAttribute('width')),count:Number(e.nextElementSibling.textContent)})));
  check(tag+' MIC cells and counts',cells.every(c=>c.width>=22)&&JSON.stringify(cells.map(c=>c.count))===JSON.stringify(micReport.mic.confusion.map(c=>c[2])),cells);
  const contrast=await heat.evaluate(region=>{
    const rgb=s=>(s.match(/[\d.]+/g)||[]).map(Number);
    const lum=c=>c.slice(0,3).map(v=>{v/=255;return v<=.04045?v/12.92:((v+.055)/1.055)**2.4}).reduce((a,v,i)=>a+v*[.2126,.7152,.0722][i],0);
    const surface=getComputedStyle(document.documentElement).getPropertyValue('--paper-pure').trim().split(' ').map(Number);
    const background=region.querySelector('svg>rect');
    const base=rgb(getComputedStyle(background).fill),alpha=base[3]??1;
    const under=base.slice(0,3).map((v,i)=>v*alpha+surface[i]*(1-alpha));
    return [...region.querySelectorAll('svg g>rect')].map(r=>{
      const fill=rgb(getComputedStyle(r).fill),fg=rgb(getComputedStyle(r.nextElementSibling).fill),a=fill[3]??1;
      const bg=fill.slice(0,3).map((v,i)=>v*a+under[i]*(1-a)),x=lum(fg),y=lum(bg);
      return {count:r.nextElementSibling.textContent,ratio:(Math.max(x,y)+.05)/(Math.min(x,y)+.05)};
    });
  });
  check(tag+' MIC text contrast',contrast.every(r=>r.ratio>=4.5),contrast);
  const probeId=manifest.panels.find(p=>p.panel==='salmonella_qrdr').dataset;
  const probeReport=JSON.parse(fs.readFileSync(path.join(FIXTURE,'datasets',probeId+'.json')));
  const probe=probeReport.panels.find(p=>p.id==='salmonella_qrdr');
  await page.goto(BASE+'/datasets?dataset='+probeId);
  const noCall=page.getByText(/The probes give no call for/);await noCall.waitFor();await noCall.scrollIntoViewIfNeeded();
  check(tag+' probe no-call details',await noCall.textContent().then(t=>t.includes(`${probe.no_call.resistant} resistant and ${probe.no_call.susceptible} susceptible`)),await noCall.textContent());
  await shot(page,tag+'-probe-no-calls');
}

async function states(browser, width, height, theme) {
  const tag=width+'-'+theme;
  const context=await browser.newContext({viewport:{width,height},reducedMotion:'reduce'});
  await context.addInitScript(t=>localStorage.setItem('grm-theme',t),theme);
  const page=await context.newPage();
  page.on('pageerror',e=>errors.push({tag,state:true,message:e.message}));
  let release;
  const gate=new Promise(resolve=>release=resolve);
  await page.route('**/data/index.json',async route=>{await gate;await route.continue();});
  await page.goto(BASE,{waitUntil:'domcontentloaded'});
  await page.getByText('Loading the report manifest',{exact:true}).waitFor();
  check(tag+' loading aria-busy',await page.locator('[aria-busy="true"]').count()===1,null);
  await shot(page,tag+'-loading');release();
  await page.getByRole('heading',{name:'Overview',exact:true}).waitFor();
  await page.unroute('**/data/index.json');
  await page.route('**/data/index.json',route=>route.fulfill({status:500,body:'unavailable'}));
  await page.reload();await page.getByRole('alert').waitFor();await shot(page,tag+'-manifest-error');
  await page.unroute('**/data/index.json');
  await page.route('**/data/datasets/escherichia_coli__ceftriaxone.json',route=>route.fulfill({status:404,body:'not found'}));
  await page.goto(BASE+'/datasets?dataset=escherichia_coli__ceftriaxone');
  await page.getByRole('alert').waitFor();await shot(page,tag+'-dataset-error');
  await page.unroute('**/data/datasets/escherichia_coli__ceftriaxone.json');
  const empty={...manifest,datasets:[],panels:[]};
  await page.route('**/data/index.json',route=>route.fulfill({json:empty}));
  await page.goto(BASE);await page.getByText('No dataset has a report yet',{exact:true}).waitFor();await shot(page,tag+'-empty-manifest');
  await page.unroute('**/data/index.json');
  const mixed=structuredClone(manifest);mixed.datasets[0].evaluation.revision=manifest.evaluation_revision-1;
  await page.route('**/data/index.json',route=>route.fulfill({json:mixed}));
  await page.goto(BASE+'/datasets');await page.getByRole('button',{name:'Provisional report'}).waitFor();
  await page.getByRole('button',{name:'Provisional report'}).click();await shot(page,tag+'-mixed-audit');
  check(tag+' mixed audit count',await page.getByText('1 of '+expected.datasets+' evaluations are provisional.',{exact:true}).count()===1,null);
  await page.unroute('**/data/index.json');
  for(const [id,word] of [['enterococcus_faecium__chloramphenicol','Fold lacks a class'],['salmonella_enterica__temocillin','Too few lineages'],['klebsiella_pneumoniae__cefpodoxime_ceftriaxone','No rule'],['acinetobacter_baumannii__doxycycline','No MIC']]) {
    const row=manifest.datasets.find(d=>d.id===id);
    await page.goto(BASE+'/datasets');await page.getByRole('searchbox',{name:'Search datasets'}).fill(row.antibiotic+' '+row.species);
    check(tag+' '+id+' state',await page.getByText(word,{exact:true}).count()>0,word);
    await shot(page,tag+'-state-'+id);
  }
  await context.close();
}


async function normalMotion(browser) {
  const context=await browser.newContext({viewport:{width:390,height:844},reducedMotion:'no-preference'});
  const page=await context.newPage();
  await page.goto(BASE);await page.getByRole('button',{name:'Present',exact:true}).click();
  const dialog=page.getByRole('dialog',{name:'GRM findings presentation'});await dialog.waitFor();
  check('normal motion media preference',await page.evaluate(()=>!matchMedia('(prefers-reduced-motion: reduce)').matches),null);
  for(let slide=0;slide<6;slide++) {
    await dialog.getByRole('combobox',{name:'Go to slide'}).selectOption(String(slide));
    if(slide===5)await dialog.getByRole('button',{name:'Show sequence'}).waitFor();
    const before=await dialog.locator('main').innerText();
    await page.evaluate(()=>Promise.allSettled(document.getAnimations().map(a=>a.finished)));
    check('normal motion slide '+(slide+1)+' final values from first paint',before===await dialog.locator('main').innerText(),null);
  }
  await shot(page,'normal-motion-rule');
  await context.close();
}

async function boundaries(browser) {
  const context=await browser.newContext({viewport:{width:900,height:900},reducedMotion:'reduce',hasTouch:true});
  const page=await context.newPage();
  for(const width of [899,900,1439,1440,1680]) {
    await page.setViewportSize({width,height:900});await page.goto(BASE+'/datasets');await page.getByRole('searchbox',{name:'Search datasets'}).waitFor();
    const docked=await page.locator('[data-filter-rail="docked"]').count();
    check('rail at '+width,docked===(width>=1440?1:0),docked);
    await shot(page,'boundary-'+width);
  }
  for(const [width,zoom] of [[720,200],[360,400]]) {
    await page.setViewportSize({width,height:900});await page.goto(BASE+'/datasets');await page.getByRole('searchbox',{name:'Search datasets'}).waitFor();
    check(zoom+'% equivalent reflow',await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),await page.evaluate(()=>({width:innerWidth,scroll:document.documentElement.scrollWidth})));
    await shot(page,'reflow-'+zoom);
  }
  await page.setViewportSize({width:390,height:844});await page.emulateMedia({forcedColors:'active'});await page.goto(BASE+'/datasets');await page.getByRole('searchbox',{name:'Search datasets'}).waitFor();
  const targets=await page.evaluate(()=>[...document.querySelectorAll('button,input[type="search"]')].filter(e=>!e.disabled&&e.getBoundingClientRect().width>0&&e.getBoundingClientRect().top<innerHeight&&e.getBoundingClientRect().bottom>0).map(e=>({name:e.getAttribute('aria-label')??e.textContent,width:e.getBoundingClientRect().width,height:e.getBoundingClientRect().height})).filter(e=>e.width<44||e.height<44));
  check('coarse pointer 44px targets',targets.length===0,targets);
  await page.getByRole('button',{name:'Random',exact:true}).focus();await shot(page,'forced-colors-focus');
  check('forced-colors status glyph',await page.locator('[data-verdict]').first().textContent().then(t=>/[✓×]/.test(t)),null);
  await context.close();
}
(async () => {
  const browser = await chromium.launch({ executablePath:'/root/.cache/ms-playwright/chromium_headless_shell-1248/chrome-headless-shell-linux64/chrome-headless-shell', headless:true, args:['--no-sandbox','--disable-dev-shm-usage'] });
  try {
    for (const [width,height] of [[1280,720],[1440,900],[390,844]]) for (const theme of ['light','dark']) {
      const context = await browser.newContext({ viewport:{width,height}, colorScheme:theme, reducedMotion:'reduce' });
      await context.addInitScript(theme => localStorage.setItem('grm-theme',theme),theme);
      const page = await context.newPage();
      const cdp = await context.newCDPSession(page);
      await cdp.send('Emulation.setDeviceMetricsOverride',{width,height,deviceScaleFactor:1,mobile:false});
      page.on('pageerror',error=>errors.push({theme,width,message:error.message}));
      page.on('response',r=>{if(r.status()>=400&&!r.url().endsWith('/data/research.json'))errors.push({theme,width,url:r.url(),status:r.status()});});
      const tag = width + '-' + theme;
      try {
        await page.goto(BASE);
        await page.getByRole('heading',{name:'Overview',exact:true}).waitFor();
        await page.getByRole('button',{name:'Present',exact:true}).click();
        const dialog=page.getByRole('dialog',{name:'GRM findings presentation'});
        await dialog.waitFor();
        for(let slide=0;slide<6;slide++) {
          await dialog.getByRole('combobox',{name:'Go to slide'}).selectOption(String(slide));
          if(slide===5)await dialog.getByRole('button',{name:'Show sequence'}).waitFor();
          await shot(page,tag+'-slide-'+(slide+1));
          const m=await page.evaluate(measure);
          check(tag+' slide '+(slide+1)+' true viewport',m.viewport[0]===width&&m.viewport[1]===height,m.viewport);
          check(tag+' slide '+(slide+1)+' fits',m.main.scrollHeight<=m.main.height+1&&m.main.scrollWidth<=m.main.width+1&&m.section.scrollHeight<=m.section.height+1,m);
          check(tag+' slide '+(slide+1)+' figure does not overlap caption',m.contentBottom<=m.caption.y-8&&m.fit.scrollHeight<=m.fit.height+1,{contentBottom:m.contentBottom,caption:m.caption,fit:m.fit});
          check(tag+' slide '+(slide+1)+' labels inside',m.textOutside.length===0,m.textOutside);
          const count=[0,expected.paired.length,expected.catalog.length,expected.panels.length,expected.mic.length,0][slide];
          check(tag+' slide '+(slide+1)+' all dataset marks',m.marks===count,{actual:m.marks,expected:count});
          await referenceChecks(page,tag,slide);
        }
        await dialog.getByRole('button',{name:'Show sequence'}).click();
        await page.getByRole('dialog',{name:'Rule sequence'}).waitFor();await shot(page,tag+'-sequence');
        await page.keyboard.press('Escape');
        check(tag+' sequence focus return',(await activeName(page))==='Show sequence',await activeName(page));
        const closeButton=dialog.getByRole('button',{name:'Close',exact:true});
        const selector=dialog.getByRole('combobox',{name:'Go to slide'});
        await closeButton.focus();await page.keyboard.press('Home');
        check(tag+' Home reaches first slide',await selector.inputValue()==='0',null);
        await page.keyboard.press('ArrowRight');
        const chart=dialog.getByRole('listbox');await chart.focus();await page.keyboard.press('ArrowRight');
        check(tag+' chart owns arrow keys',await selector.inputValue()==='1'&&await chart.getByRole('option').getAttribute('aria-posinset')==='2',null);
        await shot(page,tag+'-chart-focus');
        await page.keyboard.press('Escape');
        check(tag+' chart Escape keeps Present open',await dialog.count()===1&&await dialog.getByRole('tooltip').count()===0,null);
        await page.keyboard.press('Tab');
        check(tag+' chart allows Tab exit',!await chart.evaluate(e=>e.contains(document.activeElement)),await activeName(page));
        await closeButton.focus();await page.keyboard.press('End');
        check(tag+' End reaches last slide',await selector.inputValue()==='5',null);
        await page.keyboard.press('ArrowLeft');
        check(tag+' ArrowLeft goes back',await selector.inputValue()==='4',null);
        await page.keyboard.press('PageDown');
        check(tag+' PageDown goes forward',await selector.inputValue()==='5',null);
        await page.keyboard.press('Escape');
        check(tag+' Present focus return',(await activeName(page))==='Present',await activeName(page));
        if(width===1280){await context.close();continue;}
        for(const [route,title] of [['/','Overview'],['/datasets','Datasets'],['/species/escherichia_coli','Escherichia coli'],['/methods','Methods']]) {
          await page.goto(BASE+route);
          await page.getByRole('heading',{name:title,exact:true}).first().waitFor();
          if(route.startsWith('/species'))await page.getByRole('grid').waitFor();
          if(route==='/methods')await page.getByRole('heading',{name:'Antibiotic name resolution'}).waitFor();
          await shot(page,tag+'-'+title.split(' ')[0].toLowerCase());
          check(tag+' '+title+' page overflow',await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),await page.evaluate(()=>({width:innerWidth,scroll:document.documentElement.scrollWidth})));
        }
        await page.goto(BASE+'/datasets?dataset=escherichia_coli__ceftriaxone');
        await page.locator('#ds-learners').waitFor();
        await shot(page,tag+'-inspector');
        await interactions(page,tag,width);
        await states(browser,width,height,theme);
      } catch(error){ errors.push({theme,width,message:error.stack}); }
      await context.close();
    }
    await boundaries(browser);
    await normalMotion(browser);
  } finally { await browser.close(); }
  const receipt={source:{fixture:FIXTURE,generated:manifest.generated_at,revision:manifest.evaluation_revision,datasets:expected.datasets,panelRows:manifest.panels.length},independent:expected.stats,checks,screenshots,errors,passed:checks.every(c=>c.pass)&&errors.length===0};
  fs.writeFileSync(path.join(OUT,'receipt.json'),JSON.stringify(receipt,null,2));
  console.log(JSON.stringify({passed:receipt.passed,checks:checks.length,failures:checks.filter(c=>!c.pass).map(c=>({name:c.name,actual:c.actual})),errors,screenshots:screenshots.length},null,2));
  process.exitCode=receipt.passed?0:1;
})().catch(e=>{console.error(e);process.exitCode=1;});
