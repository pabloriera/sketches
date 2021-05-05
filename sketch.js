const particles = [];
regions = null;
const region_size = 100;
const width_ = 1200;
const height_ = 800;
const vel_init = 0.0;
const cyclic = true;
const n_particles = 651;
const draw_numbers = false;
const init_grid = true;
var vfrom = null;
var vto = null;

settings = 
{
  'sigma': 10.0,
  'epsilon': 0.005,
  'dt': 0.25,
  'maxforce': 0.5,
  'r0': 0.0,
  'drag': 0.0,
  'regions': false,
  'lines': false,
  'radious': 8,
  'line alpha':0.04,
  'influence':1.0
}

function setup() {
  createCanvas(width_, height_);

  regions = new Regions(width_/region_size,height_/region_size);
  
  for (let i = 0; i < n_particles; i++) {
    let p = new Particle(i);
    particles.push(p);  particles[i].update();
  }
  
  regions.setup(particles);
  
  frameRate(1000);
  
  colorMode(RGB);
  vfrom = color(180, 155, 32);
  vto = color(255, 61, 139);
  
  let gui = new dat.GUI();
  gui.add(settings, 'influence', 0.01, 1.0, 0.05);
  gui.add(settings, 'sigma', 1, 30.0, 0.05);
  gui.add(settings, 'r0', 0.0, 20.0, 0.1);
  gui.add(settings, 'epsilon', 0.001, 0.5, 0.001);
  gui.add(settings, 'maxforce', 0.01, 5.0);  
  gui.add(settings, 'drag', 0.0, 10.0, 0.001);
  gui.add(settings, 'dt', 0.01, 1.0);
  gui.addFolder('Vis');
  gui.add(settings, 'regions');
  gui.add(settings, 'lines');
  gui.add(settings, 'radious', 1, 20.0);
  gui.add(settings, 'line alpha', 0.04, 1.0,0.01);
  gui.add({ 'grid':function(){ grid() }},'grid');
  gui.add({ 'shake':function(){ shake() }},'shake');
  gui.add({ 'kick':function(){ kick() }},'kick');
  
}

function draw() {

  for (let i = particles.length - 1; i >= 0; i--) 
    particles[i].update(1);
    
  for (let i = particles.length - 1; i >= 0; i--){
    particles[i].update(2);
    if (particles[i].finished()) {
      // remove this particle
      // particles.splice(i, 1);
    }
  }
  
//   settings.drag = map(mouseX,0,width_,0,10,true);
//   settings.epsilon = map(mouseY,0,height_,0.07,0,true);
  
  // if(frameCount%3==0)
  {
    background(0);

    for (let i = particles.length - 1; i >= 0; i--)
      particles[i].show(); 

    if(settings.regions){
      stroke(128);  
      for (let i =0; i<width_; i+=regions.width) {
        line(i,0,i,height_)    
      }
      for (let i =0; i<height_; i+=regions.height) {
        line(0,i,width_,i)
      }
    }
  }
}

function grid(){
    for (let i = 0; i < n_particles; i++) {
      particles[i].setup();
    }
}

function shake(){
    for (let i = 0; i < n_particles; i++) {
      particles[i].vx+=randomGaussian(0,2);
      particles[i].vy+=randomGaussian(0,2);
    }
}

function kick(){
  
  for (let j =0; j<int(random(10)) ;j++){
    let i = int(random(n_particles));
    let ph = random(0,2*PI);
    let v = randomGaussian(50,20);
    particles[i].vx+=v*Math.cos(ph);
    particles[i].vy+=v*Math.sin(ph);
  }
}


class Regions{
  
  constructor(n,m){
    this.total = n*m;
    this.m = m;
    this.n = n;
    this.width = width_/n;
    this.height = height_/m;
    this.region = new p5.TypedDict();
    for (let i = 0; i<this.total; i++) {
     this.region.set(i,[]);
    }
  }
  
  get_neighbors(r){ 
    var nei = this.region.get(r);
    return nei;
  }
  
  xytoix(i,j){   
    if (cyclic)
    {
        i = i % this.n
        j = j % this.m
        if (i<0)
          i+=this.n;
        if (j<0)
          j+=this.m;
        return i+j*this.n
    }
    else
    {
      if(i<0 || j<0 || i>=this.n || j>=this.m )
        return -1
      else
        return i+j*this.n
    }
        
  }
  
  get_all_neighbors(r){ 
    var out = []
    let rx = r % this.n;
    let ry = int(r/this.n);
    
    for (let i in [0,1,2])
      for (let j in [0,1,2])
      {
        i = int(i)
        j = int(j)
        // console.log(r,i-1,j-1,rx,ry)
        let rs = this.xytoix(rx+i-1,ry+j-1)
        // console.log(r,rx+i-1,ry+j-1,i-1,j-1,rx,ry,rs)
        // if (rs>=this.total)
          // console.log(r,rx+i-1,ry+j-1,i-1,j-1,rx,ry,rs)
        // console.log(r,rs,rx,ry,i-1,j-1);

        if (rs>=0)
        {
          var nei = this.region.get(rs);
          // console.log(nei);
          out = out.concat(nei)
        }
      }

    return out;
  }
  
  get_region(x,y){
    return int(x/this.width)+int(y/this.height)*this.n;
  }
  
  set_neighbors(i,r){
    var nei = this.region.get(r);
    nei.push(i);
  }
  
  remove_neighbors(i,r){
    var nei = this.region.get(r);
    let ni = nei.indexOf(i);
    nei.splice(ni,1);    
  }
  
  setup(particles){
    for (let i = 0; i<particles.length; i++) {
      let r = this.get_region(particles[i].x,particles[i].y);
      this.set_neighbors(i,r);
      particles[i].region = r;
    }
  }
}

function lennard_jones(r)
{
  return -48 * settings.epsilon * ( Math.pow(settings.sigma / r, 13) - 0.5 * Math.pow(settings.sigma/r, 7) )
}

class Particle {

  constructor(i) {
    this.i = i;
    this.setup()
  }
  setup(){
    if(init_grid)
    {
      let a = Math.sqrt(n_particles/(width_/height_));
      this.x = ((this.i % 31)+0.5)*width_/31;
      this.y = (int(this.i/31)+0.5)*height_/21;
    }
    else
    {
      this.x = random(0, width_);
      this.y = random(0, height_);
    }
    this.region = regions.get_region(this.x,this.y);
    this.vel = random(0,vel_init);
    this.angle = random(0,2*PI);
    this.vx = this.vel*sin(this.angle);
    this.vy = this.vel*cos(this.angle);
    this.ax = 0;
    this.ay = 0;
    this.alpha = 255;
    this.neighbors = [];
  }

  finished() {
    return this.alpha < 0;
  }
  
  leapfrog(part)
  {
    if(part==1){
      this.vx = this.vx + 0.5*settings.dt*this.ax;
      this.vy = this.vy + 0.5*settings.dt*this.ay;
      this.x = this.x + settings.dt*this.vx;
      this.y = this.y + settings.dt*this.vy;
    }
    else if(part==2)
    {
      this.vx = this.vx + 0.5*settings.dt*this.ax;
      this.vy = this.vy + 0.5*settings.dt*this.ay;      
    }
  }

  forces(){    
    let nei = this.neighbors;
    this.ax = 0;
    this.ay = 0;
    let forcex = 0;
    let forcey = 0;
    for (let i = 0; i<nei.length; i++) {
      let force = 0;
      let dx = particles[nei[i]].x-this.x;
      let dy = particles[nei[i]].y-this.y;
      let rr = dx*dx+dy*dy;
      let r = Math.sqrt(rr);
      if (r < settings.influence*region_size)
      {
        if (r>0){
          force = lennard_jones(r-settings.r0);    
          force = max(min(force,settings.maxforce),-settings.maxforce);
          forcex += force*dx;
          forcey += force*dy;
        }
        // else
        // {
        //   let v = createVector(this.vx,this.vy);
        //   let f = v.normalize().mult(-1).mult(settings.maxforce);
        //   forcex = f.x
        //   forcey = f.y;
        //   console.log('pasa');
        // }          
      }
      
      
    }
    this.ax += forcex;
    this.ay += forcey;

    let v = createVector(this.vx,this.vy);
    let f = v.normalize().mult(-1).mult(Math.pow(v.mag(),2)).mult(settings.drag);
    forcex = f.x;
    forcey = f.y;
    
    this.ax += forcex;
    this.ay += forcey;    
  }
    
  update(part) {    
    if(part==1)
      this.leapfrog(part);
    else if(part==2)
    {
      this.forces();
      this.leapfrog(part);

      if (cyclic){
        this.x = this.x % width_;
        if (this.x<0) this.x+=width_;
        this.y = this.y % height_;
        if (this.y<0) this.y+=height_;  
      }
      else
      {
        if (this.x<0 || this.y<0 || this.x>width_ || this.y>height_)
          this.alpha = -1;
      }

      if (this.alpha>0)
      {
        let rnow = regions.get_region(this.x,this.y);
        if (rnow!=this.region)
        {
          regions.remove_neighbors(this.i, this.region);
          regions.set_neighbors(this.i, rnow);
          this.region = rnow;
        }    
      }

      this.neighbors = regions.get_all_neighbors(this.region);
      let index = this.neighbors.indexOf(this.i);
      if (index > -1) {
         this.neighbors.splice(index, 1);
      }
    }
  }

  show() {
    if (!this.finished())
    {      
      if(settings.lines){
        stroke(255,255,255, settings['line alpha']*255);
        let nei = this.neighbors; 
        for (let i = 0; i<nei.length; i++) {
          let dx = particles[nei[i]].x-this.x;
          let dy = particles[nei[i]].y-this.y;
          let r = Math.sqrt(dx*dx+dy*dy)
          if(r<settings.influence*region_size)
          {
            let overx = 0;
            let overy = 0;
            if(dx>region_size)
              overx = -width_;
            else if(dx<-region_size)
              overx = width_;

            if(dy>region_size)
              overy = -height_;
            else if(dy<-region_size)
              overy = height_;

            line(this.x,this.y, 
                 particles[nei[i]].x+overx, 
                 particles[nei[i]].y+overy)
          }
        }
      }

      
      let vel = Math.sqrt(this.vx*this.vx+this.vy*this.vy)/2;
      stroke(lerpColor(vfrom, vto, vel/2));
      fill(lerpColor(vfrom, vto, vel), this.alpha);
      ellipse(this.x, this.y, settings.radious);
      if (draw_numbers){
        fill(0,255);
        textSize(10)
        text(this.i, this.x*0.999, this.y*1.01);}
      }
    
  }

}