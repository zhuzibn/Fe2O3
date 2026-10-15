/* Notebook-based drawing only. These are NOT fractional crystallographic
 * coordinates, physical lattice constants, or an input to the MATLAB model.
 * The author explicitly chose an illustrative view rather than a sourced CIF.
 */
'use strict';
const el = id => document.getElementById(id);
const canvas = el('crystal'), ctx = canvas.getContext('2d');
const H = 3.5, SQRT3 = Math.sqrt(3);
const t1 = [1.5,SQRT3/2,0], t2 = [0,SQRT3,0];
// Illustrative separation of the four sites and three column sequences.
const heights = [0.06,0.20,0.56,0.70], phases = [0,1/3,2/3];
// Notebook labels run A–D from top to bottom in each column sequence.
const columnLevels = phases.map(phase=>heights.map(z=>1-((z+phase)%1)).sort((a,b)=>b-a));
// One outline per illustrated layer: pair the nearby staggered marker levels.
// When both perimeter groups occur, its vertices follow their respective heights.
const markerLevels = columnLevels.flatMap((levels,group)=>levels.map(z=>({group,z}))).sort((a,b)=>b.z-a.z);
const guideLayers = Array.from({length:markerLevels.length/2},(_,i)=>{
  const perimeter=markerLevels.slice(2*i,2*i+2).filter(level=>level.group!==2);
  return [0,1].map(group=>(perimeter.find(level=>level.group===group)||perimeter[0]).z);
});
const names = ['A','B','C','D'], groupColors = ['#283746','#167fb5','#cd791d'];
let atoms = [], edges = [], layerEdges = [], picked = null, hitTargets = [], bounds, framePending = false;
let pointer = null;
const key = p => p.map(x=>Math.round(x*1e7)).join(',');
const add = (a,b) => a.map((v,i)=>v+b[i]);
function repeats() { return ['repeat-a','repeat-b','repeat-c'].map(id=>Number(el(id).value)); }
function build() {
  const [na,nb,nc] = repeats(), atomMap = new Map(), edgeMap = new Map(), layerEdgeMap = new Map();
  const vertex = Array.from({length:6},(_,v)=>[Math.cos(v*Math.PI/3),Math.sin(v*Math.PI/3),0]);
  for(let a=0;a<na;a++) for(let b=0;b<nb;b++) for(let c=0;c<nc;c++) {
    const origin = [a*t1[0]+b*t2[0],a*t1[1]+b*t2[1],c*H];
    const corners=vertex.map(v=>add(origin,v));
    const edge=(p,q)=>{
      const id=[key(p),key(q)].sort().join('|');
      if(!edgeMap.has(id)) edgeMap.set(id,[p,q]);
    };
    for(let v=0;v<6;v++) {
      const p=corners[v],q=corners[(v+1)%6];
      edge(p,q); edge(add(p,[0,0,H]),add(q,[0,0,H])); edge(p,add(p,[0,0,H]));
      // A single layer loop joins its staggered perimeter positions.
      for(const layer of guideLayers) {
        const lp=add(p,[0,0,H*layer[v%2]]),lq=add(q,[0,0,H*layer[(v+1)%2]]);
        const id=[key(lp),key(lq)].sort().join('|');
        if(!layerEdgeMap.has(id)) layerEdgeMap.set(id,[lp,lq]);
      }
    }
    const columns=[...vertex.map((p,v)=>({p,group:v%2})),{p:[0,0,0],group:2}];
    for(let v=0;v<columns.length;v++) {
      const column=columns[v];
      for(let n=0;n<4;n++) {
        const z=columnLevels[column.group][n]*H;
        const pos=add(origin,[column.p[0],column.p[1],z]);
        const id=key(pos), owner=[a,b,c];
        if(atomMap.has(id)) { atomMap.get(id).owners.push(owner); continue; }
        atomMap.set(id,{id,pos,group:column.group,letter:names[n],sign:n===0||n===3?1:-1,owners:[owner],label:a===0&&b===0&&c===0&&(v===0||v===1||v===6)});
      }
    }
  }
  atoms=[...atomMap.values()]; edges=[...edgeMap.values()]; layerEdges=[...layerEdgeMap.values()];
  const points=edges.flat();
  const low=[0,1,2].map(i=>Math.min(...points.map(p=>p[i])));
  const high=[0,1,2].map(i=>Math.max(...points.map(p=>p[i])));
  bounds={center:low.map((v,i)=>(v+high[i])/2),span:high.map((v,i)=>v-low[i])};
  picked=null;
  el('selection').textContent='Select an Fe sphere to inspect its illustration label.';
  el('counts').textContent=`${na*nb*nc} drawn prism${na*nb*nc===1?'':'s'} · ${atoms.length} Fe markers, including boundary columns.`;
  schedule();
}
function view() {
  const yaw=Number(el('yaw').value)*Math.PI/180, tilt=Number(el('tilt').value)*Math.PI/180;
  const width=canvas.clientWidth,height=canvas.clientHeight;
  const radius=Math.hypot(...bounds.span)/2;
  const scale=Math.min(width-90,height-100)/(2*radius)*Number(el('zoom').value)/100;
  const rotate=p=>{
    const [x,y,z]=p;
    const x1=x*Math.cos(yaw)-y*Math.sin(yaw), y1=x*Math.sin(yaw)+y*Math.cos(yaw);
    return [x1,-z*Math.cos(tilt)+y1*Math.sin(tilt),y1*Math.cos(tilt)+z*Math.sin(tilt)];
  };
  const project=p=>{
    const q=rotate(p.map((v,i)=>v-bounds.center[i]));
    return [width/2+q[0]*scale,height/2+q[1]*scale,q[2]];
  };
  return {width,height,scale,project,rotate};
}
function line(p,q,color,width=1) {
  ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.strokeStyle=color;ctx.lineWidth=width;ctx.stroke();
}
function arrow(p,q,color,width=2) {
  line(p,q,color,width);
  const angle=Math.atan2(q[1]-p[1],q[0]-p[0]);
  const size=width===2?7:5;
  ctx.beginPath();ctx.moveTo(q[0],q[1]);
  ctx.lineTo(q[0]-size*Math.cos(angle-.45),q[1]-size*Math.sin(angle-.45));
  ctx.lineTo(q[0]-size*Math.cos(angle+.45),q[1]-size*Math.sin(angle+.45));
  ctx.closePath();ctx.fillStyle=color;ctx.fill();
}
function draw() {
  framePending=false;
  const v=view(), dpr=Math.min(window.devicePixelRatio||1,2);
  const width=Math.round(v.width*dpr),height=Math.round(v.height*dpr);
  if(canvas.width!==width||canvas.height!==height) {canvas.width=width;canvas.height=height;}
  ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,v.width,v.height);
  const objects=[];
  if(el('frames').checked) for(const [p,q] of edges) {
    const pp=v.project(p),qq=v.project(q);
    objects.push({depth:(pp[2]+qq[2])/2,paint:()=>line(pp,qq,'#b9b4a9',1)});
  }
  if(el('frames').checked) for(const [p,q] of layerEdges) {
    const pp=v.project(p),qq=v.project(q);
    objects.push({depth:(pp[2]+qq[2])/2,paint:()=>{
      ctx.save();ctx.setLineDash([4,4]);line(pp,qq,'#b9b4a9',1);ctx.restore();
    }});
  }
  hitTargets=[];
  for(const atom of atoms) {
    const p=v.project(atom.pos), r=Math.max(2.4,Math.min(11,v.scale*.10));
    const color=el('color-mode').value==='columns'?groupColors[atom.group]:'#b76a25';
    objects.push({depth:p[2],paint:()=>{
      if(el('spins').checked) {
        const end=v.project(add(atom.pos,[atom.sign*.40,0,0]));
        arrow(p,end,atom.sign===1?'#b51b87':'#21854b',atoms.length>500?1.4:2);
      }
      const gradient=ctx.createRadialGradient(p[0]-r*.35,p[1]-r*.4,r*.07,p[0],p[1],r);
      gradient.addColorStop(0,'#fff6cb');gradient.addColorStop(.32,color);gradient.addColorStop(1,'#35291e');
      ctx.beginPath();ctx.arc(p[0],p[1],r,0,Math.PI*2);ctx.fillStyle=gradient;ctx.fill();
      if(atom.id===picked) {ctx.beginPath();ctx.arc(p[0],p[1],r+4,0,Math.PI*2);ctx.strokeStyle='#ca920e';ctx.lineWidth=3;ctx.stroke();}
      if(el('labels').checked&&atom.label) {
        ctx.font='13px Georgia';ctx.fillStyle='#243746';
        ctx.fillText(atom.letter+['₁','₂','₃'][atom.group],p[0]+r+7,p[1]-5);
      }
      hitTargets.push({atom,x:p[0],y:p[1],r:Math.max(r+4,8)});
    }});
  }
  // Painter's order for the orthographic illustration, not a solid-surface renderer.
  objects.sort((a,b)=>a.depth-b.depth).forEach(o=>o.paint());
  const anchor=[58,v.height-52];
  for(const [name,vec] of [['a₁',t1],['a₂',t2],['c',[0,0,1.8]]]) {
    const q=v.rotate(vec),end=[anchor[0]+q[0]*24,anchor[1]+q[1]*24];
    arrow(anchor,end,'#445a69',1.4);ctx.font='13px Georgia';ctx.fillStyle='#243746';ctx.fillText(name,end[0]+4,end[1]-4);
  }
  ctx.fillStyle='#657078';ctx.font='12px "Segoe UI", sans-serif';
  ctx.fillText('NOTEBOOK SCHEMATIC · NOT TO SCALE',18,26);
}
function schedule() {if(!framePending) {framePending=true;requestAnimationFrame(draw);}}
function camera(yaw,tilt) {el('yaw').value=yaw;el('tilt').value=tilt;el('zoom').value=100;schedule();}
function inspect(x,y) {
  const hit=[...hitTargets].reverse().find(p=>Math.hypot(p.x-x,p.y-y)<=p.r);
  picked=hit?hit.atom.id:null;
  el('selection').textContent=hit?`${hit.atom.letter}${hit.atom.group+1} · illustration column group ${hit.atom.group+1} · drawn prism indices ${hit.atom.owners.map(o=>'('+o.map(n=>n+1).join(',')+')').join(', ')}. No MATLAB site index assigned.`:'Select an Fe sphere to inspect its illustration label.';
  schedule();
}
canvas.addEventListener('pointerdown',event=>{
  if(event.button!==0) return;
  canvas.focus({preventScroll:true});canvas.setPointerCapture(event.pointerId);
  pointer={id:event.pointerId,x:event.clientX,y:event.clientY,startX:event.clientX,startY:event.clientY,moved:false};
});
canvas.addEventListener('pointermove',event=>{
  if(!pointer||pointer.id!==event.pointerId) return;
  const dx=event.clientX-pointer.x,dy=event.clientY-pointer.y;
  pointer.moved ||= Math.hypot(event.clientX-pointer.startX,event.clientY-pointer.startY)>4;
  el('yaw').value=((Number(el('yaw').value)+dx*.5+540)%360)-180;
  el('tilt').value=Math.max(-85,Math.min(85,Number(el('tilt').value)+dy*.4));
  pointer.x=event.clientX;pointer.y=event.clientY;schedule();
});
canvas.addEventListener('pointerup',event=>{
  if(!pointer||pointer.id!==event.pointerId) return;
  if(!pointer.moved) {const box=canvas.getBoundingClientRect();inspect(event.clientX-box.left,event.clientY-box.top);}
  pointer=null;
});
canvas.addEventListener('pointercancel',()=>{pointer=null;});
canvas.addEventListener('keydown',event=>{
  if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(event.key)) return;
  event.preventDefault();
  if(event.key==='ArrowLeft'||event.key==='ArrowRight') el('yaw').value=((Number(el('yaw').value)+(event.key==='ArrowLeft'?-5:5)+540)%360)-180;
  else el('tilt').value=Math.max(-85,Math.min(85,Number(el('tilt').value)+(event.key==='ArrowUp'?5:-5)));
  schedule();
});
for(const id of ['repeat-a','repeat-b','repeat-c']) el(id).addEventListener('change',build);
for(const id of ['yaw','tilt','zoom','frames','spins','labels','color-mode']) el(id).addEventListener('input',schedule);
el('single').addEventListener('click',()=>{for(const id of ['repeat-a','repeat-b','repeat-c']) el(id).value='1';build();});
el('larger').addEventListener('click',()=>{for(const id of ['repeat-a','repeat-b','repeat-c']) el(id).value='2';build();});
el('reset').addEventListener('click',()=>camera(-25,15));
el('top').addEventListener('click',()=>camera(-25,85));
el('side').addEventListener('click',()=>camera(0,0));
new ResizeObserver(schedule).observe(canvas);
build();
