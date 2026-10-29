/* Documentation snapshot of notebook pp. 8, 10–11 and current MATLAB indexing.
 * These arrays describe implemented offsets, not physical crystal coordinates.
 * Keep the table synchronized with structure-exchange-guide.md and the sources.
 */
'use strict';
const palette = {green:'#3f7c35', orange:'#d56b16', blue:'#1378b8', purple:'#7b3ea1', red:'#d92d27', black:'#171717'};
const letters = {green:'G',orange:'O',blue:'B',purple:'P',red:'R',black:'K'};
const P = [[0,1],[1,-1],[-1,-1]], M = [[0,-1],[1,1],[-1,1]];
const P2 = [[0,2],[1,0],[-1,0]], M2 = [[0,-2],[1,0],[-1,0]], O = [[0,0]];
const group = (color, offsets, layer) => offsets.map(([r,c]) => ({color,offset:[r,c,layer]}));
const sites = {
  green: {layer:1,mask:'atomtype_layer1gr',col:9,
    bonds:[group('red',O,-1),group('orange',P2,0),group('purple',M,1),[...group('blue',P2,1),...group('black',M,-1)]],
    lines:[[12,14],[18,23],[38,43],[78,86]],surface:[[7,10],[18,23],[38,43],[63,69]]},
  orange: {layer:1,mask:'atomtype_layer1or',col:11,
    bonds:[group('blue',O,1),group('green',M2,0),group('black',P,-1),[...group('purple',P,1),...group('red',M2,-1)]],
    lines:[[12,14],[25,30],[50,55],[88,96]],surface:[[7,10],[25,30],[45,48],[71,76]]},
  blue: {layer:2,mask:'atomtype_layer2blue',col:11,
    bonds:[group('orange',O,-1),group('purple',P,0),group('red',M2,1),[...group('black',P,1),...group('green',M2,-1)]],
    lines:[[104,106],[110,115],[131,136],[151,159]]},
  purple: {layer:2,mask:'atomtype_layer2p',col:12,
    bonds:[group('black',O,1),group('blue',M,0),group('green',P,-1),[...group('red',P,1),...group('orange',M,-1)]],
    lines:[[104,106],[118,123],[138,143],[161,169]]},
  red: {layer:3,mask:'atomtype_layer3red',col:9,
    bonds:[group('green',O,1),group('black',M,0),group('blue',P2,-1),[...group('purple',M,-1),...group('orange',P2,1)]],
    lines:[[182,184],[195,200],[208,213],[258,266]],surface:[[177,180],[195,200],[208,213],[241,246]]},
  black: {layer:3,mask:'atomtype_layer3black',col:12,
    bonds:[group('purple',O,-1),group('red',P,0),group('orange',M,1),[...group('blue',M,-1),...group('green',P,1)]],
    lines:[[182,184],[188,193],[220,225],[248,256]],surface:[[177,180],[188,193],[215,218],[234,239]]}
};
// dmi.m reuses J3/J4 sums. Signs depend on the actual layer, not motif color.
// Weights are in joules; these labels do not reconstruct physical bond vectors.
const dmiDefaults = {3:1.6e-23,4:1.28e-23};
function dmiBonds(color,k) {
  const layerSign=1-2*(k%2);
  return [3,4].flatMap(shell=>sites[color].bonds[shell-1].map(b=>{
    const coefficient=layerSign*(shell===3?1:-1);
    return {...b,shell,coefficient,weight:coefficient*dmiDefaults[shell]};
  }));
}
const signed = n => n>0?'+':'−';
const weightLabel = b => `${signed(b.coefficient)}${Math.abs(b.weight).toExponential().replace('e-23',' × 10⁻²³')} J`;
const $ = id => document.getElementById(id);
const esc = value => String(value).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
let selected = -1;
function colorAt(r,c,k) {
  if (r<1 || r>9 || c<1 || c>17 || k<1 || k>6) return null;
  const motif=(k-1)%3, remainder=((c-(r%2 ? 1 : 3))%4+4)%4;
  if (motif===0) return remainder===0?'green':remainder===2?'orange':null;
  if (motif===1) return remainder===2?'blue':remainder===3?'purple':null;
  return remainder===0?'red':remainder===3?'black':null;
}
function shiftName([r,c,k]) {
  const suffix = `${r>0?'u':r<0?'d':''}${r===0?'':Math.abs(r)}${c>0?'l':c<0?'r':''}${c===0?'':Math.abs(c)}`;
  const variable = r===0 && c===0 ? 'mmxtmp' : `mmx_${suffix}`;
  return `${variable}(:,:,i${k===0?'':k>0?`+${k}`:k})`;
}
function state() {
  const color=$('target').value, site=sites[color], k=Number($('boundary').value), surface=k===1||k===6;
  const edge=$('edge').value==='left';
  const r=edge && !['green','red'].includes(color)?4:5;
  const c=edge?(['purple','black'].includes(color)?2:1):site.col;
  const dmi=$('shell').value==='DMI', shell=dmi?null:Number($('shell').value.slice(1))-1;
  const bonds=dmi?dmiBonds(color,k):site.bonds[shell];
  const neighbors=bonds.map((b,index)=>{
    const pos=b.offset.map((n,j)=>n+[r,c,k][j]);
    const missingLayer=pos[2]<1 || pos[2]>6;
    const missingEdge=pos[0]<1 || pos[0]>9 || pos[1]<1 || pos[1]>17;
    return {...b,index,pos,active:!missingLayer&&!missingEdge,reason:missingLayer?'Outside stack':missingEdge?'Empty padding':'Included',shift:shiftName(b.offset)};
  });
  return {color,site,surface,r,c,k,shell,dmi,neighbors};
}
function circle(x,y,color,label,active=true,highlight=false,target=false) {
  return `${target?`<circle cx="${x}" cy="${y}" r="16" fill="none" stroke="#bd8516" stroke-width="3"/>`:''}<circle cx="${x}" cy="${y}" r="${highlight?13:10}" fill="${active?palette[color]:'#fffdf8'}" stroke="${active?palette[color]:'#a34b35'}" stroke-width="${highlight?4:2}" ${active?'':'stroke-dasharray="3 3"'}/><text x="${x}" y="${y+4}" text-anchor="middle" style="fill:${active?'white':'#843d2a'};font-size:11px">${label}</text>`;
}
function drawBonds(s) {
  const project=(r,c,k)=>[260+c*55+r*24,195+k*115+r*20];
  let svg='';
  for (let dk=-1;dk<=1;dk++) {
    const y=195+dk*115, layer=s.k+dk, absent=layer<1||layer>6;
    svg+=`<path d="M55 ${y-20} L390 ${y-20} L470 ${y+30} L135 ${y+30} Z" fill="${absent?'#eeeae3':'#f3f0e8'}" stroke="#cfc6b5" ${absent?'stroke-dasharray="5 5"':''}/><text x="14" y="${y-30}" style="font-size:13px">k=${layer}${absent?' · absent layer':''}${s.dmi&&dk!==0?(dk<0?' · physical +z':' · physical −z'):''}</text>`;
  }
  for (const b of s.neighbors) {
    const [x,y]=project(...b.offset);
    const bondColor=s.dmi?(b.shell===3?'#125d88':'#ae6633'):palette[b.color];
    svg+=`<line x1="260" y1="195" x2="${x}" y2="${y}" stroke="${b.index===selected?'#a46700':b.active?bondColor:'#b7ada0'}" stroke-width="${b.index===selected?4:1.5}" ${b.active?'':'stroke-dasharray="4 4"'}/>`;
  }
  for (const b of s.neighbors) {
    const [x,y]=project(...b.offset);
    svg+=circle(x,y,b.color,b.index+1,b.active,b.index===selected);
  }
  if(s.dmi) for(const b of s.neighbors.filter(b=>b.active)) {
    const [x,y]=project(...b.offset);
    svg+=`<text class="dmi-coefficient" x="${x}" y="${y+25}" text-anchor="middle" style="font-size:12px;paint-order:stroke;stroke:#fffdf8;stroke-width:3px">${signed(b.coefficient)}D${b.shell}</text>`;
  }
  svg+=circle(260,195,s.color,letters[s.color],true,false,true);
  svg+='<text x="14" y="375" style="font-size:12px">r increases diagonally ↘ · c increases → · k increases downward</text>';
  $('bond-view').innerHTML=svg;
}
function drawArrays(s) {
  let html='';
  for(let dk=-1;dk<=1;dk++) {
    const layer=s.k+dk, absent=layer<1||layer>6;
    let svg='';
    for(let dc=-3;dc<=3;dc++) {
      if(s.c+dc < -1 || s.c+dc > 19) continue; // Two padded columns only.
      svg+=`<text x="${48+(dc+3)*32}" y="15" text-anchor="middle" style="font-size:11px">${s.c+dc}</text>`;
    }
    for(let dr=-2;dr<=2;dr++) {
      const r=s.r+dr,y=38+(dr+2)*32;
      svg+=`<text x="12" y="${y+4}" style="font-size:11px">${r}</text>`;
      for(let dc=-3;dc<=3;dc++) {
        const c=s.c+dc,x=48+(dc+3)*32,padding=r<1||r>9||c<1||c>17;
        if(r < -1 || r > 11 || c < -1 || c > 19) continue; // Outside allocated array.
        if(padding && $('padding').value==='hide') continue;
        const b=s.neighbors.find(n=>n.pos[0]===r&&n.pos[1]===c&&n.pos[2]===layer);
        const target=dr===0&&dc===0&&dk===0;
        const color=colorAt(r,c,layer);
        svg+=`<rect x="${x-16}" y="${y-16}" width="32" height="32" fill="${padding||absent?'#e5e1da':b&&b.index===selected?'#f4e4ba':'white'}" stroke="#cfc6b5"/>`;
        if(target) svg+=circle(x,y,s.color,letters[s.color],true,false,true);
        else if(b) svg+=circle(x,y,b.color,b.index+1,b.active,b.index===selected);
        else if(color) svg+=`<text x="${x}" y="${y+5}" text-anchor="middle" style="fill:${palette[color]};font-size:16px;font-weight:bold">${letters[color]}</text>`;
        else svg+=`<text x="${x}" y="${y+5}" text-anchor="middle" style="fill:#858178">${absent?'–':padding?'×':'·'}</text>`;
      }
    }
    html+=`<div class="array-layer"><h3>k=${layer}${absent?' · absent':''}</h3><svg viewBox="0 0 275 185" role="img" aria-label="Array window at layer ${layer}">${svg}</svg></div>`;
  }
  $('arrays').innerHTML=html;
}
function sourceExcerpt(source,start,end,highlight='') {
  return source.split('\n').slice(start-1,end).map((line,j)=>{
    let body=esc(line.replace(/^\s{12}/,''));
    if(highlight) body=body.replaceAll(esc(highlight),`<mark>${esc(highlight)}</mark>`);
    return `<span class="ln">${start+j}  </span>${body}`;
  }).join('\n');
}
function render() {
  const s=state(), count=s.neighbors.filter(b=>b.active).length;
  if(s.dmi&&selected<0) selected=s.neighbors.findIndex(b=>b.active);
  $('summary').textContent=`${s.color[0].toUpperCase()+s.color.slice(1)} (${s.r},${s.c},${s.k}) · ${s.dmi?'DMI · D3 + D4':`J${s.shell+1}`} · ${count} of ${s.neighbors.length} listed neighbors present`;
  $('dmi-note').hidden=!s.dmi;
  if(s.dmi) {
    const layerSign=1-2*(s.k%2);
    $('dmi-layer').textContent=`k=${s.k}: (−1)^k = ${signed(layerSign)}1. J3 weight: ${weightLabel(s.neighbors[0])}; J4 weight: ${weightLabel(s.neighbors[3])}. Present: ${s.neighbors.filter(b=>b.active&&b.shell===3).length}/3 J3 and ${s.neighbors.filter(b=>b.active&&b.shell===4).length}/6 J4 neighbors.`;
  }
  $('mask-note').textContent=`Target selector: ${s.site.mask}. Padded target: (${s.r+2},${s.c+2},${s.k}). ${s.site.layer===2?'Blue/purple layers do not terminate the illustrated six-layer stack.':''}`;
  drawBonds(s); drawArrays(s);
  $('neighbors').innerHTML=s.neighbors.map(b=>`<tr class="${selected===b.index?'selected':''}"><td><button type="button" data-neighbor="${b.index}" aria-pressed="${selected===b.index}" aria-label="Highlight neighbor ${b.index+1}">${b.index+1}</button></td><td><span class="swatch" style="background:${palette[b.color]}"></span>${b.color}</td><td>${s.dmi?`J${b.shell} · ${b.active?`${signed(b.coefficient)}D${b.shell}<br>${weightLabel(b)}`:'omitted'}`:`J${s.shell+1}`}</td><td>(${b.pos.join(', ')})</td><td>(${b.offset.join(', ')})</td><td><code>${esc(b.shift)}</code></td><td>${b.reason}</td></tr>`).join('');
  const selectedBond=s.neighbors[selected];
  const ranges=s.surface?s.site.surface:s.site.lines;
  const [start,end]=s.dmi?[1,10]:ranges[s.shell];
  const filename=s.dmi?'dmi.m':'exchangej_.m';
  $('source-location').innerHTML=`<a href="../${filename}">${filename}</a> · lines ${start}–${end} · ${s.dmi?'shell weighting':s.surface?'surface branch':'interior branch'}`;
  $('source-code').innerHTML=sourceExcerpt(s.dmi?dmiSource:exchangeSource,start,end,selectedBond?(s.dmi?`mmxtmpJ${selectedBond.shell}`:selectedBond.shift):'');
  $('dmi-neighbor-source').hidden=!s.dmi;
  if(s.dmi) {
    const shells=selectedBond?[selectedBond.shell]:[3,4];
    $('neighbor-source').innerHTML=shells.map(shell=>{
      const [a,b]=ranges[shell-1];
      return `<p class="caption"><a href="../exchangej_.m">exchangej_.m</a> · J${shell} · lines ${a}–${b} · ${s.surface?'surface':'interior'} branch</p><pre tabindex="0">${sourceExcerpt(exchangeSource,a,b,selectedBond?.shift)}</pre>`;
    }).join('');
  }
  $('explanation').textContent=s.dmi?`The ${s.color} mask supplies J3/J4 neighbor sums in exchangej_.m. dmi.m weights them by (−1)^k D3 and −(−1)^k D4. Each weight W contributes +W·sy/μs to Bx and −W·sx/μs to By; Bz is zero. Select a neighbor to trace its shift and shell sum in both source excerpts. Dashed neighbors contribute zero because their layer or lateral site is absent.`:s.surface&&count===0?'The selected surface branch omits this target contribution because its neighboring layer is absent. The dashed neighbors show the interior list for comparison.':`The color mask selects ${s.color} targets. Shifted arrays gather the neighboring spins; their sum is the J${s.shell+1} neighbor sum before multiplication by the exchange constant. ${s.color==='blue'&&s.shell===3?'The y/z temporary names say “black”, but their mask and use correctly select blue targets.':''}`;
  const page=s.dmi?8:s.site.layer===3?11:10;
  $('page-label').textContent=`page ${page}`;
  $('pdf-link').href=`../Alpha-Fe2O3.pdf#page=${page}`;
  $('notebook-image').src=`Alpha-Fe2O3-pages/page-${String(page).padStart(2,'0')}.jpg`;
  $('notebook-image').alt=s.dmi?'Historical notebook page 8: J3/J4 note; current code uses those neighbor sums with explicit shell weights and layer signs':`Notebook page ${page}: original structure sketches, square-array encoding, and exchange annotations`;
}
function changeSelection(event) {
  const motif=sites[$('target').value].layer, previous=Number($('boundary').value);
  const layers=[motif,motif+3];
  $('boundary').innerHTML=layers.map(k=>`<option value="${k}">k=${k} · ${k===1||k===6?'surface':'interior'} · ${k%2?'odd':'even'}</option>`).join('');
  $('boundary').value=layers.includes(previous)?previous:motif===1?4:motif;
  if(!event||['target','shell'].includes(event.target.id)) selected=-1;
  render();
}
function initialize() {
  $('dmi-parameters').textContent=dmiParameterSource;
  $('dmi-field').textContent=dmiFieldSource;
  $('legend').innerHTML=Object.keys(palette).map(color=>`<span><i class="dot" style="background:${palette[color]}"></i>${letters[color]} · ${color}</span>`).join('');
  for(const id of ['target','shell','boundary','edge','padding']) $(id).addEventListener('change',changeSelection);
  $('neighbors').addEventListener('click',event=>{
    const button=event.target.closest('button[data-neighbor]');
    if(!button) return;
    selected=Number(button.dataset.neighbor);
    render();
    $('neighbors').querySelector(`button[data-neighbor="${selected}"]`).focus({preventScroll:true});
  });
  changeSelection();
}
// Verbatim source snapshot. Line numbers above refer to this text.
const exchangeSource = "\nfor i=1:natomH\n    switch mod(i-1,3)\n        case 0\n            %% calculate J1\n            % green atom in the 1st layer don't have J1\n            if i==1\n                mmxtmpJ1(:,:,i)=atomtype_layer1or.*mmxtmp(:,:,i+1);\n                mmytmpJ1(:,:,i)=atomtype_layer1or.*mmytmp(:,:,i+1);\n                mmztmpJ1(:,:,i)=atomtype_layer1or.*mmztmp(:,:,i+1);\n            else\n                mmxtmpJ1(:,:,i)=atomtype_layer1gr.*mmxtmp(:,:,i-1)+atomtype_layer1or.*mmxtmp(:,:,i+1);\n                mmytmpJ1(:,:,i)=atomtype_layer1gr.*mmytmp(:,:,i-1)+atomtype_layer1or.*mmytmp(:,:,i+1);\n                mmztmpJ1(:,:,i)=atomtype_layer1gr.*mmztmp(:,:,i-1)+atomtype_layer1or.*mmztmp(:,:,i+1);\n            end\n            %% calculate J2\n            %J2 experienced by green atoms\n            shift_matx_gr=atomtype_layer1gr.*(mmx_l2(:,:,i)+mmx_d1(:,:,i)...\n                +mmx_u1(:,:,i));\n            shift_maty_gr=atomtype_layer1gr.*(mmy_l2(:,:,i)+mmy_d1(:,:,i)...\n                +mmy_u1(:,:,i));\n            shift_matz_gr=atomtype_layer1gr.*(mmz_l2(:,:,i)+mmz_d1(:,:,i)...\n                +mmz_u1(:,:,i));\n            %J2 experienced by orange atoms\n            shift_matx_or=atomtype_layer1or.*(mmx_r2(:,:,i)+mmx_d1(:,:,i)+...\n                mmx_u1(:,:,i));\n            shift_maty_or=atomtype_layer1or.*(mmy_r2(:,:,i)+mmy_d1(:,:,i)+...\n                mmy_u1(:,:,i));\n            shift_matz_or=atomtype_layer1or.*(mmz_r2(:,:,i)+mmz_d1(:,:,i)+...\n                mmz_u1(:,:,i));\n\n            mmxtmpJ2(:,:,i)=shift_matx_gr+shift_matx_or;\n            mmytmpJ2(:,:,i)=shift_maty_gr+shift_maty_or;\n            mmztmpJ2(:,:,i)=shift_matz_gr+shift_matz_or;\n\n            %% calculate J3\n            %J3 experienced by green atoms\n            shift_matx_gr=atomtype_layer1gr.*(mmx_r1(:,:,i+1)+mmx_u1l1(:,:,i+1)+...\n                mmx_d1l1(:,:,i+1));\n            shift_maty_gr=atomtype_layer1gr.*(mmy_r1(:,:,i+1)+mmy_u1l1(:,:,i+1)+...\n                mmy_d1l1(:,:,i+1));\n            shift_matz_gr=atomtype_layer1gr.*(mmz_r1(:,:,i+1)+mmz_u1l1(:,:,i+1)+...\n                mmz_d1l1(:,:,i+1));\n            %J3 experienced by orange atoms\n            if i==1\n                shift_matx_or=0;\n                shift_maty_or=0;\n                shift_matz_or=0;\n            else\n                shift_matx_or=(atomtype_layer1or.*(mmx_l1(:,:,i-1)+mmx_u1r1(:,:,i-1)+...\n                    mmx_d1r1(:,:,i-1)));\n                shift_maty_or=(atomtype_layer1or.*(mmy_l1(:,:,i-1)+mmy_u1r1(:,:,i-1)+...\n                    mmy_d1r1(:,:,i-1)));\n                shift_matz_or=(atomtype_layer1or.*(mmz_l1(:,:,i-1)+mmz_u1r1(:,:,i-1)+...\n                    mmz_d1r1(:,:,i-1)));\n            end\n            mmxtmpJ3(:,:,i)=shift_matx_gr+shift_matx_or;\n            mmytmpJ3(:,:,i)=shift_maty_gr+shift_maty_or;\n            mmztmpJ3(:,:,i)=shift_matz_gr+shift_matz_or;\n\n            %% calculate J4\n            %J4 experienced by green atoms\n            if i==1\n                shift_matx_gr=atomtype_layer1gr.*(mmx_l2(:,:,i+1)+...\n                    mmx_d1(:,:,i+1)+mmx_u1(:,:,i+1));\n                shift_maty_gr=atomtype_layer1gr.*(mmy_l2(:,:,i+1)+...\n                    mmy_d1(:,:,i+1)+mmy_u1(:,:,i+1));\n                shift_matz_gr=atomtype_layer1gr.*(mmz_l2(:,:,i+1)+...\n                    mmz_d1(:,:,i+1)+mmz_u1(:,:,i+1));\n                %J4 experienced by orange atoms\n                shift_matx_or=atomtype_layer1or.*(mmx_l1(:,:,i+1)+...\n                    mmx_u1r1(:,:,i+1)+mmx_d1r1(:,:,i+1));\n                shift_maty_or=atomtype_layer1or.*(mmy_l1(:,:,i+1)+...\n                    mmy_u1r1(:,:,i+1)+mmy_d1r1(:,:,i+1));\n                shift_matz_or=atomtype_layer1or.*(mmz_l1(:,:,i+1)+...\n                    mmz_u1r1(:,:,i+1)+mmz_d1r1(:,:,i+1));\n            else\n                shift_matx_gr=atomtype_layer1gr.*(mmx_l2(:,:,i+1)+...\n                    mmx_d1(:,:,i+1)+mmx_u1(:,:,i+1)+(mmx_r1(:,:,i-1)+...\n                    mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1)));\n                shift_maty_gr=atomtype_layer1gr.*(mmy_l2(:,:,i+1)+...\n                    mmy_d1(:,:,i+1)+mmy_u1(:,:,i+1)+(mmy_r1(:,:,i-1)+...\n                    mmy_u1l1(:,:,i-1)+mmy_d1l1(:,:,i-1)));\n                shift_matz_gr=atomtype_layer1gr.*(mmz_l2(:,:,i+1)+...\n                    mmz_d1(:,:,i+1)+mmz_u1(:,:,i+1)+(mmz_r1(:,:,i-1)+...\n                    mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1)));\n                %J4 experienced by orange atoms\n                shift_matx_or=atomtype_layer1or.*(mmx_l1(:,:,i+1)+...\n                    mmx_u1r1(:,:,i+1)+mmx_d1r1(:,:,i+1)+(mmx_r2(:,:,i-1)+...\n                    mmx_u1(:,:,i-1)+mmx_d1(:,:,i-1)));\n                shift_maty_or=atomtype_layer1or.*(mmy_l1(:,:,i+1)+...\n                    mmy_u1r1(:,:,i+1)+mmy_d1r1(:,:,i+1)+ (mmy_r2(:,:,i-1)+...\n                    mmy_u1(:,:,i-1)+mmy_d1(:,:,i-1)));\n                shift_matz_or=atomtype_layer1or.*(mmz_l1(:,:,i+1)+...\n                    mmz_u1r1(:,:,i+1)+mmz_d1r1(:,:,i+1)+(mmz_r2(:,:,i-1)+...\n                    mmz_u1(:,:,i-1)+mmz_d1(:,:,i-1)));\n            end\n            mmxtmpJ4(:,:,i)=shift_matx_gr+shift_matx_or;\n            mmytmpJ4(:,:,i)=shift_maty_gr+shift_maty_or;\n            mmztmpJ4(:,:,i)=shift_matz_gr+shift_matz_or;\n\n        case 1\n            %% calculate J1\n            mmxtmpJ1(:,:,i)=atomtype_layer2p.*mmxtmp(:,:,i+1)+atomtype_layer2blue.*mmxtmp(:,:,i-1);\n            mmytmpJ1(:,:,i)=atomtype_layer2p.*mmytmp(:,:,i+1)+atomtype_layer2blue.*mmytmp(:,:,i-1);\n            mmztmpJ1(:,:,i)=atomtype_layer2p.*mmztmp(:,:,i+1)+atomtype_layer2blue.*mmztmp(:,:,i-1);\n\n            %% calculate J2\n            %J2 experienced by blue atoms\n            shift_matx_blue=atomtype_layer2blue.*(mmx_l1(:,:,i)+mmx_u1r1(:,:,i)+...\n                mmx_d1r1(:,:,i));\n            shift_maty_blue=atomtype_layer2blue.*(mmy_l1(:,:,i)+mmy_u1r1(:,:,i)+...\n                mmy_d1r1(:,:,i));\n            shift_matz_blue=atomtype_layer2blue.*(mmz_l1(:,:,i)+mmz_u1r1(:,:,i)+...\n                mmz_d1r1(:,:,i));\n\n            %J2 experienced by purple atoms\n            shift_matx_purple=atomtype_layer2p.*(mmx_r1(:,:,i)+mmx_u1l1(:,:,i)+...\n                mmx_d1l1(:,:,i));\n            shift_maty_purple=atomtype_layer2p.*(mmy_r1(:,:,i)+mmy_u1l1(:,:,i)+...\n                mmy_d1l1(:,:,i));\n            shift_matz_purple=atomtype_layer2p.*(mmz_r1(:,:,i)+mmz_u1l1(:,:,i)+...\n                mmz_d1l1(:,:,i));\n\n            mmxtmpJ2(:,:,i)=shift_matx_blue+shift_matx_purple;\n            mmytmpJ2(:,:,i)=shift_maty_blue+shift_maty_purple;\n            mmztmpJ2(:,:,i)=shift_matz_blue+shift_matz_purple;\n\n            %% calculate J3\n            %J3 experienced by blue atoms\n            shift_matx_blue=atomtype_layer2blue.*(mmx_r2(:,:,i+1)+mmx_u1(:,:,i+1)+...\n                mmx_d1(:,:,i+1));\n            shift_maty_blue=atomtype_layer2blue.*(mmy_r2(:,:,i+1)+mmy_u1(:,:,i+1)+...\n                mmy_d1(:,:,i+1));\n            shift_matz_blue=atomtype_layer2blue.*(mmz_r2(:,:,i+1)+mmz_u1(:,:,i+1)+...\n                mmz_d1(:,:,i+1));\n            %J3 experienced by purple atoms\n            shift_matx_purple=atomtype_layer2p.*(mmx_l1(:,:,i-1)+mmx_u1r1(:,:,i-1)+...\n                mmx_d1r1(:,:,i-1));\n            shift_maty_purple=atomtype_layer2p.*(mmy_l1(:,:,i-1)+mmy_u1r1(:,:,i-1)+...\n                mmy_d1r1(:,:,i-1));\n            shift_matz_purple=atomtype_layer2p.*(mmz_l1(:,:,i-1)+mmz_u1r1(:,:,i-1)+...\n                mmz_d1r1(:,:,i-1));\n\n            mmxtmpJ3(:,:,i)=shift_matx_blue+shift_matx_purple;\n            mmytmpJ3(:,:,i)=shift_maty_blue+shift_maty_purple;\n            mmztmpJ3(:,:,i)=shift_matz_blue+shift_matz_purple;\n\n            %% calculate J4\n            %J4 experienced by blue atoms\n            shift_matx_blue=atomtype_layer2blue.*(mmx_l1(:,:,i+1)+...\n                mmx_u1r1(:,:,i+1)+mmx_d1r1(:,:,i+1)+mmx_r2(:,:,i-1)+...\n                mmx_u1(:,:,i-1)+mmx_d1(:,:,i-1));\n            shift_maty_black=atomtype_layer2blue.*(mmy_l1(:,:,i+1)+...\n                mmy_u1r1(:,:,i+1)+mmy_d1r1(:,:,i+1)+mmy_r2(:,:,i-1)+...\n                mmy_u1(:,:,i-1)+mmy_d1(:,:,i-1));\n            shift_matz_black=atomtype_layer2blue.*(mmz_l1(:,:,i+1)+...\n                mmz_u1r1(:,:,i+1)+mmz_d1r1(:,:,i+1)+mmz_r2(:,:,i-1)+...\n                mmz_u1(:,:,i-1)+mmz_d1(:,:,i-1));\n            %J4 experienced by purple atoms\n            shift_matx_purple=atomtype_layer2p.*(mmx_l1(:,:,i+1)+...\n                mmx_u1r1(:,:,i+1)+mmx_d1r1(:,:,i+1)+mmx_r1(:,:,i-1)+...\n                mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1));\n            shift_maty_purple=atomtype_layer2p.*(mmy_l1(:,:,i+1)+...\n                mmy_u1r1(:,:,i+1)+mmy_d1r1(:,:,i+1)+mmy_r1(:,:,i-1)+...\n                mmy_u1l1(:,:,i-1)+mmy_d1l1(:,:,i-1));\n            shift_matz_purple=atomtype_layer2p.*(mmz_l1(:,:,i+1)+...\n                mmz_u1r1(:,:,i+1)+mmz_d1r1(:,:,i+1)+mmz_r1(:,:,i-1)+...\n                mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1));\n\n            mmxtmpJ4(:,:,i)=shift_matx_blue+shift_matx_purple;\n            mmytmpJ4(:,:,i)=shift_maty_black+shift_maty_purple;\n            mmztmpJ4(:,:,i)=shift_matz_black+shift_matz_purple;\n\n        case 2\n            %% calculate J1\n            if i==natomH\n                mmxtmpJ1(:,:,i)=atomtype_layer3black.*mmxtmp(:,:,i-1);\n                mmytmpJ1(:,:,i)=atomtype_layer3black.*mmytmp(:,:,i-1);\n                mmztmpJ1(:,:,i)=atomtype_layer3black.*mmztmp(:,:,i-1);\n            else\n                mmxtmpJ1(:,:,i)=(atomtype_layer3red.*mmxtmp(:,:,i+1))+atomtype_layer3black.*mmxtmp(:,:,i-1);\n                mmytmpJ1(:,:,i)=(atomtype_layer3red.*mmytmp(:,:,i+1))+atomtype_layer3black.*mmytmp(:,:,i-1);\n                mmztmpJ1(:,:,i)=(atomtype_layer3red.*mmztmp(:,:,i+1))+atomtype_layer3black.*mmztmp(:,:,i-1);\n            end\n            %% calculate J2\n            %J2 experienced by black atoms\n            shift_matx_black=atomtype_layer3black.*(mmx_l1(:,:,i)+mmx_u1r1(:,:,i)+...\n                mmx_d1r1(:,:,i));\n            shift_maty_black=atomtype_layer3black.*(mmy_l1(:,:,i)+mmy_u1r1(:,:,i)+...\n                mmy_d1r1(:,:,i));\n            shift_matz_black=atomtype_layer3black.*(mmz_l1(:,:,i)+mmz_u1r1(:,:,i)+...\n                mmz_d1r1(:,:,i));\n            %J2 experienced by red atoms\n            shift_matx_red=atomtype_layer3red.*(mmx_r1(:,:,i)+mmx_u1l1(:,:,i)+...\n                mmx_d1l1(:,:,i));\n            shift_maty_red=atomtype_layer3red.*(mmy_r1(:,:,i)+mmy_u1l1(:,:,i)+...\n                mmy_d1l1(:,:,i));\n            shift_matz_red=atomtype_layer3red.*(mmz_r1(:,:,i)+mmz_u1l1(:,:,i)+...\n                mmz_d1l1(:,:,i));\n\n            mmxtmpJ2(:,:,i)=shift_matx_black+shift_matx_red;\n            mmytmpJ2(:,:,i)=shift_maty_black+shift_maty_red;\n            mmztmpJ2(:,:,i)=shift_matz_black+shift_matz_red;\n\n            %% calculate J3\n            %J3 for red atoms\n            shift_matx_red=atomtype_layer3red.*(mmx_l2(:,:,i-1)+mmx_u1(:,:,i-1)+...\n                mmx_d1(:,:,i-1));\n            shift_maty_red=atomtype_layer3red.*(mmy_l2(:,:,i-1)+mmy_u1(:,:,i-1)+...\n                mmy_d1(:,:,i-1));\n            shift_matz_red=atomtype_layer3red.*(mmz_l2(:,:,i-1)+mmz_u1(:,:,i-1)+...\n                mmz_d1(:,:,i-1));\n            %J3 for black atoms\n            if i==natomH\n                shift_matx_black=0;\n                shift_maty_black=0;\n                shift_matz_black=0;\n            else\n                shift_matx_black=atomtype_layer3black.*(mmx_r1(:,:,i+1)+mmx_u1l1(:,:,i+1)+...\n                    mmx_d1l1(:,:,i+1));\n                shift_maty_black=atomtype_layer3black.*(mmy_r1(:,:,i+1)+mmy_u1l1(:,:,i+1)+...\n                    mmy_d1l1(:,:,i+1));\n                shift_matz_black=atomtype_layer3black.*(mmz_r1(:,:,i+1)+mmz_u1l1(:,:,i+1)+...\n                    mmz_d1l1(:,:,i+1));\n            end\n            mmxtmpJ3(:,:,i)=shift_matx_red+shift_matx_black;\n            mmytmpJ3(:,:,i)=shift_maty_red+shift_maty_black;\n            mmztmpJ3(:,:,i)=shift_matz_red+shift_matz_black;\n\n            %% calculate J4\n            %J4 experienced by black atoms\n            if i==natomH\n                shift_matx_black=atomtype_layer3black.*(mmx_r1(:,:,i-1)+...\n                    mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1));\n                shift_maty_black=atomtype_layer3black.*(mmy_r1(:,:,i-1)+...\n                    mmy_u1l1(:,:,i-1)+mmy_d1l1(:,:,i-1));\n                shift_matz_black=atomtype_layer3black.*(mmz_r1(:,:,i-1)+...\n                    mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1));\n                %J4 experienced by red atoms\n                shift_matx_red=atomtype_layer3red.*(mmx_r1(:,:,i-1)+...\n                    mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1));\n                shift_maty_red=atomtype_layer3red.*(mmy_r1(:,:,i-1)+mmy_u1l1(:,:,i-1)+...\n                    mmy_d1l1(:,:,i-1));\n                shift_matz_red=atomtype_layer3red.*(mmz_r1(:,:,i-1)+...\n                    mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1));\n            else\n                shift_matx_black=atomtype_layer3black.*((mmx_l1(:,:,i+1)+...\n                    mmx_u1r1(:,:,i+1)+mmx_d1r1(:,:,i+1))+mmx_r1(:,:,i-1)+...\n                    mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1));\n                shift_maty_black=atomtype_layer3black.*((mmy_l1(:,:,i+1)+...\n                    mmy_u1r1(:,:,i+1)+mmy_d1r1(:,:,i+1))+mmy_r1(:,:,i-1)+...\n                    mmy_u1l1(:,:,i-1)+mmy_d1l1(:,:,i-1));\n                shift_matz_black=atomtype_layer3black.*((mmz_l1(:,:,i+1)+...\n                    mmz_u1r1(:,:,i+1)+mmz_d1r1(:,:,i+1))+mmz_r1(:,:,i-1)+...\n                    mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1));\n                %J4 experienced by red atoms\n                shift_matx_red=atomtype_layer3red.*((mmx_d1(:,:,i+1)+...\n                    mmx_u1(:,:,i+1)+mmx_l2(:,:,i+1))+mmx_r1(:,:,i-1)+...\n                    mmx_u1l1(:,:,i-1)+mmx_d1l1(:,:,i-1));\n                shift_maty_red=atomtype_layer3red.*((mmy_d1(:,:,i+1)+...\n                    mmy_u1(:,:,i+1)+mmy_l2(:,:,i+1))+mmy_r1(:,:,i-1)+mmy_u1l1(:,:,i-1)+...\n                    mmy_d1l1(:,:,i-1));\n                shift_matz_red=atomtype_layer3red.*((mmz_d1(:,:,i+1)+...\n                    mmz_u1(:,:,i+1)+ mmz_l2(:,:,i+1))+mmz_r1(:,:,i-1)+...\n                    mmz_u1l1(:,:,i-1)+mmz_d1l1(:,:,i-1));\n            end\n            mmxtmpJ4(:,:,i)=shift_matx_black+shift_matx_red;\n            mmytmpJ4(:,:,i)=shift_maty_black+shift_maty_red;\n            mmztmpJ4(:,:,i)=shift_matz_black+shift_matz_red;\n        otherwise\n            error('unexpected index');\n    end\nend\n\n";
const dmiSource = "%% Shell-resolved z-DMI neighbor sums [J].\n% Phys. Rev. B 107, 184426(2023)\n% DMI of J1, J3 and J4\nfor i=1:natomH\n    dmi_layer_sign=1-2*mod(i,2);\n    dmi_sumx(:,:,i)=dmi_layer_sign.*(Dsim1.*mmxtmpJ1(:,:,i)+Dsim3.*mmxtmpJ3(:,:,i)-Dsim4.*mmxtmpJ4(:,:,i));\n    dmi_sumy(:,:,i)=dmi_layer_sign.*(Dsim1.*mmytmpJ1(:,:,i)+Dsim3.*mmytmpJ3(:,:,i)-Dsim4.*mmytmpJ4(:,:,i));\n    dmi_sumz(:,:,i)=dmi_layer_sign.*(Dsim1.*mmztmpJ1(:,:,i)+Dsim3.*mmztmpJ3(:,:,i)-Dsim4.*mmztmpJ4(:,:,i));\nend\nclear dmi_layer_sign\n";
const dmiParameterSource = "% Shell-resolved z-DMI [J]; J1 is retained at zero, J2 is omitted.\n% D3/D4=1.25 is a Fig.3(b) estimate; absolute values are user-selected.\nif DMIenable\n    Dsim1=0;\n    Dsim3=1.6e-23;\n    Dsim4=1.28e-23;\nelse\n    Dsim1=0;\n    Dsim3=0;\n    Dsim4=0;\nend";
const dmiFieldSource = "%% DMI\ndmi()\ndmi_sumx=dmi_sumx.*atomtype_;\ndmi_sumy=dmi_sumy.*atomtype_;\ndmi_sumz=dmi_sumz.*atomtype_;\n\nhdmi_x=dmi_sumy./muigpu;\nhdmi_y=-dmi_sumx./muigpu;\nhdmi_z=zeros(size(hex_x),'like',hex_x);";
initialize();
