import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class sketch_180727a_clusteredNetwork extends PApplet {



Network nn = new Network();
float[] testErrorList=new float[0];
float[] activeCells=new float[0];
float minCells = 100;
float err = 0;
int timer = 0;
PFont font;
//------------------------------
public void setup() {
  // load data:
  //buildData();
  importData();

  // set up NeuralNet
  nn.cellsPerCluster = 4; // set size of clusters
  nn.batchSize = 0; // set batchlearning size
  nn.squaredError = false;
  nn.memory = 1;  // lookback if the dataset is time-based!
  //nn.addLayer(4, 4);  // in - out
  nn.addLayer(4,12);
  nn.addLayer(12,6);
    nn.addLayer(6,6);
    nn.addLayer(6,6);
  nn.addLayer(6,4);

  nn.lr = 1.0e-1f; // learningRate


  size(1920, 990);
  strokeCap(RECT);
  frameRate(Integer.MAX_VALUE);
//  font = loadFont("CourierNewPSMT-48.vlw");
//  textFont(font,40);
//frameRate(1000);
  timer = millis();
  background(0);
}


//------------------------------

public void draw() {
  err = nn.learn(trainingIn, trainingOut);
  frame.setTitle("error: "+nf(100*err, 3, 2)+"% // "+nf(round(frameRate), 5, 0)+" cps lr: "+nf(nn.lr, 1, 9));

  // only update frame every x steps
  if (frameCount%100!=0)return;
  //if (frameCount%2000!=0)importData(); // shuffle input ar
  testErrorList = append(testErrorList, err);
  float[] ret;
  float cellsActive = 0, maxCells = 0;
  String output = "";
  background(0);
  pushMatrix();


  translate(300, height/3);
  for (int i=0; i<validIn.length; i++) {


    ret = nn.forward(validIn[i]);
    stroke(255);
    /*
    for (int j=0; j<inputData[i].length; j++) {
    drawLine(0, inputData[i][j]);
    }*/


    translate(0, height/3);
    //translate(5, 0);

    if(validOut[i][0]>0) stroke(200,100,0);
    if(validOut[i][2]>0) stroke(0,100,200);
    if(validOut[i][1]>0) stroke(100,0,200);
    if(validOut[i][3]>0) stroke(0,200,0);
    drawLine(0, 1);

    int index = highestIndex(ret);
    if(index==0) stroke(200,100,0);
    if(index==2) stroke(0,100,200);
    if(index==1) stroke(100,0,200);
    if(index==3) stroke(0,200,0);
    drawLine(0, ret[index]);

    translate(0, -height/3);


    for (int l=0; l<nn.layers.length; l++) {
      for (int j=0; j<nn.layers[l].clusterManager.Sy.length; j++) {
        cellsActive += nn.layers[l].clusterManager.Sy[j]? 1: 0;
        maxCells++;
        for(int k=0;k<nn.layers[l].clusterManager.clusterSize;k++)
          output += nn.layers[l].clusterManager.Sy[j]? "\u25a0":"\u25a1";
      }
      output+= "     ";
    }

    output+= "\n";
  }

  popMatrix();

  textSize(5);
  textLeading(4);
  fill(255);
  stroke(255);
  text(output, 50, 50);




  cellsActive/=validIn.length;
  maxCells/=validIn.length;
  minCells = cellsActive<minCells?cellsActive:minCells;

  stroke(255);
  displayError(testErrorList);

  fill(255);
  stroke(225);
  textSize(12);
  text("cells utilized: "+nf((float)cellsActive/maxCells*100,3,2)+" %    ("+nf((float)cellsActive,2,2)+" / "+(int)maxCells+")",10,20);
}
//------------------------------
public void drawLine(float y, float x) {
  strokeWeight(3);
  line(0, y, 0, map(x, 2, -2, -height/3, height/3));
  translate(4, 0);
}
//------------------------------
public int highestIndex(float ar[]){
  int highestIndex = -1;
  float val = -Float.MAX_VALUE;

  for(int i=0;i<ar.length;i++){
    if(ar[i] >= val){
      val = ar[i];
      highestIndex = i;
    }
  }
  return highestIndex;
}
//------------------------------
public void displayError(float[] list) {
  stroke(0, 200, 0);
  stroke(150);
  strokeWeight(1);
  for (int i=0; i<list.length-1; i++)
  line(
    map(i, 0, list.length-1, 0, width),
    map(list[i], 0, max(list), height, 0),
    map(i+1, 0, list.length-1, 0, width),
    map(list[i+1], 0, max(list), height, 0)
    );
  }
  //------------------------------
  public void keyPressed() {
    if (key=='+') nn.lr*=2;
    if (key=='-') nn.lr/=2;
  }
/*
......[INPUTS]
.......I    I
..[SWTCH]   I
..I     I   I
..I    [CLUSTER]
..I     I
...[OUT]

*/
/*
class ClusterManager {
  Cluster[] clusters;
  int clusterSize;
  boolean[] Sy;
  //------------------------------
  // setup
  ClusterManager(int in, int out, int clusterSize) {
    this.clusterSize=clusterSize;
    Sy = new boolean[out];
    int numClusters = (out/clusterSize>=1?out/clusterSize:1 ) +1;
    clusters = new Cluster[numClusters];
    for (int i=0; i<clusters.length; i++)
    clusters[i]=new Cluster(in);
    //  clusters[0].alwaysOn=true; // shouldn't need that!
  }
  //------------------------------
  void forward(float[] in) {
    int j;
    boolean thisState;
    // set Sy, by going through all the clusters
    float sumIn = sum(in);
    for (int i=0; i<clusters.length; i++) {
      thisState = clusters[i].forward(sumIn);
      j=0;
      while (j<clusterSize && i*clusterSize+j<Sy.length) {
        Sy[i*clusterSize+j] = thisState;
        j++;
      }
    }
  }
  //------------------------------
  void backward(float[] dy) {
    float ddy;
    for (int i=0; i<clusters.length; i++){
      ddy = 0;
      for(int j=0;j<clusterSize;j++){
        ddy+=dy[i/clusterSize+j];		// hmmmmm.....
      }
      //ddy/=clusterSize;
      clusters[i].backward(ddy); // I think that's not right....
    }
  }

  //------------------------------
  void update(float lr) {
    for (int i=0; i<clusters.length; i++)
    clusters[i].update(lr);
  }
  //------------------------------
  boolean isActive(int index){
    return Sy[index];
  }
  //------------------------------
  float sum(float ar[]){
    float ret = 0;
    for(int i=0;i<ar.length;i++)
    ret+=ar[i];
    return ret;
  }
  //------------------------------
}


//------------------------------
//
//------------------------------
class Cluster {
  float S0;
  float Bx, dBx;
  float Sx;
  float[] Wx,dWx;
  float W,dW;
  float threshold=1f;
	float dT;
  boolean currentState;
  boolean alwaysOn=false;
  int lr;
  //------------------------------
  Cluster(int in) {

    Wx = new float[in];
    dWx = new float[in];
    W = random(-1,1)*1f;
    for(int i=0;i<Wx.length;i++)
    Wx[i]=random(-1,1)*1;
  }
  //------------------------------
  boolean forward(float in) {
    if(alwaysOn){return true;}
    S0 = in;

    Sx = in*W+Bx;
    // add non-linearity?
    //  Sx = 2 * Sx / (1 + abs(Sx)); // softSign

    currentState = Sx > threshold ? true:false;

    return currentState;
  }
  //------------------------------
  void backward(float dy) {
    if(alwaysOn)return;

    //  dy = (2.0 / sq(1.0 + abs(Sx))) * dy; //softsign

    //if(currentState && dy<0 || !currentState && dy>0){
    if(true){

      dBx+= dy;
    //  dW+=dy;
      dW+=S0*dy;

    }
    // lr -> 1f/dy?


  }
  //------------------------------
  float sum(float ar[]){
    float ret = 0;
    for(int i=0;i<ar.length;i++)
    ret+=ar[i];
    return ret;
  }
  //------------------------------
  void update(float lr) {
    if(alwaysOn)return;

    lr*=1e-1; // adjust cluster slower than Cells!


  W+=dW*lr;
  dW=0;


  Bx += dBx*lr;
  dBx=0;

}
//------------------------------
}
*/
float[][] inputData, targetData;
float[][] trainingIn,trainingOut,validIn,validOut;

public void buildData() {
  inputData=new float[4][2];
  targetData = new float[4][1];
  // XOR
  inputData[0][0] = 0;
  inputData[0][1] = 0;
  inputData[1][0] = 0;
  inputData[1][1] = 1;
  inputData[2][0] = 1;
  inputData[2][1] = 0;
  inputData[3][0] = 1;
  inputData[3][1] = 1;
  targetData[0][0] = 1;
  targetData[1][0] = 0;
  targetData[2][0] = 0;
  targetData[3][0] = 1;

  if(true)return;

  int numSamples = 100;
  inputData=new float[numSamples][1];
  targetData = new float[numSamples][1];
  for (int i=1; i<inputData.length; i++) {
    //for (int j=0; j<inputData[i].length; j++) {
    inputData[i][0] = sin(i*.1f);
    // inputData[i][1] = cos(i*.1);
    //}
  }
  for (int i=0; i<inputData.length; i++) {
    targetData[i][0] =cos(i*.2f);
  }

  for (int i=1; i<inputData.length; i++) {
    for (int j=0; j<inputData[i].length; j++) {
      inputData[i][0] =sin(i*.1f);
      // inputData[i][1] = cos(i*.1);
    }
  }
  for (int i=0; i<inputData.length; i++) {

    targetData[i][0] = abs(sin(i*.1f)*cos(i*.1f));
  }
















  numSamples = round(15*TWO_PI);
  inputData=new float[numSamples][2];
  targetData = new float[numSamples][1];
  for (int i=1; i<inputData.length; i++) {
    inputData[i][0] = sin(radians(i*10));
    inputData[i][0] += sin(radians(i*20));
    inputData[i][1] = cos(radians(i*5));

    inputData[i][0]/=2;
  }
  for (int i=0; i<inputData.length; i++) {
    targetData[i][0] =cos(radians(i*10));
    targetData[i][0] +=sin(radians(i*5));
    targetData[i][0]/=2;
  }
}




//------------------------------------------
public void importData() {
  println("importing iris-dataset");
  float[][] f, t;
  String lines[];
  String tmp[][];

  // build validation data
  lines = loadStrings("data/iris_extended.data");
  //shuffleArray(lines);
  tmp   = new String[lines.length][0];
  for (int i=0; i<lines.length; i++)
  tmp[i] = split(lines[i], ",");



  f = new float[tmp.length][4];
  t = new float[tmp.length][4];

  for (int i=0; i<tmp.length; i++) {
    f[i][0] = PApplet.parseInt(tmp[i][0]);
    f[i][1] = PApplet.parseInt(tmp[i][1]);
    f[i][2] = PApplet.parseInt(tmp[i][2]);
    f[i][3] = PApplet.parseInt(tmp[i][3]);
    if(tmp[i][4].equals("Iris-setosa"))
    t[i][0] = 1;
    else if(tmp[i][4].equals("Iris-virginica"))
    t[i][2] = 1;
    else if(tmp[i][4].equals("Iris-mythica"))
    t[i][1] = 1;
    else if(tmp[i][4].equals("Iris-versicolor"))
    t[i][3] = 1;


  }
  trainingIn = f;
  trainingOut = t;
  validIn = f;
  validOut = t;

  /*maybe make a test set with x times the validation set but shuffled every time
  so i can maybe use lookback on the memory-cells? could be woth a shot...*/

  if(true)return;


  // build validation data
  lines = loadStrings("data/iris.data");

  tmp   = new String[lines.length][0];
  for (int i=0; i<lines.length; i++)
  tmp[i] = split(lines[i], ",");

  f = new float[tmp.length][4];
  t = new float[tmp.length][1];

  for (int i=0; i<tmp.length; i++) {
    f[i][0] = PApplet.parseInt(tmp[i][0]);
    f[i][1] = PApplet.parseInt(tmp[i][1]);
    f[i][2] = PApplet.parseInt(tmp[i][2]);
    f[i][3] = PApplet.parseInt(tmp[i][3]);
    if(tmp[i][4].equals("Iris-setosa"))
    t[i][0] = -1;
    else if(tmp[i][4].equals("Iris-versicolor"))
    t[i][0] = 0;
    else if(tmp[i][4].equals("Iris-virginica"))
    t[i][0] = 1;
  }
  validIn = f;
  validOut = t;

}

public void shuffleArray(float[][] ar1, float[][] ar2)
{
  if(ar1.length!=ar2.length)return;

  println("shuffling");
  int index; float[] temp1,temp2;
  java.util.Random random = new java.util.Random();
  for (int i = ar1.length - 1; i > 0; i--)
  {
    index = random.nextInt(i + 1);
    temp1 = ar1[index];
    temp2 = ar2[index];
    ar1[index] = ar1[i];
    ar1[i] = temp1;

    ar2[index] = ar2[i];
    ar2[i] = temp2;
  }
}


public void shuffleArray(String[] array)
{
  println("shuffling");
  int index; String temp;
  Random random = new Random();
  for (int i = array.length - 1; i > 0; i--)
  {
    index = random.nextInt(i + 1);
    temp = array[index];
    array[index] = array[i];
    array[i] = temp;
  }
}
//------------------------------
class Layer {
  MemoryCell cells[];
  float[] Sy;
  float[] dSx;
  int inputs, outputs;
  ClusterManager clusterManager;
  int layerLevel;
  //------------------------------
  // init layer
  Layer(int _in, int _out, int clusterSize,int memory,int level) {
    inputs=_in;
    outputs = _out;
    layerLevel = level;
    Sy = new float[outputs];
    dSx = new float[inputs];
    cells = new MemoryCell[outputs];

    for (int i=0; i<outputs; i++) {
      cells[i] = new MemoryCell(inputs,memory,level);
    }
    clusterManager=new ClusterManager(inputs,outputs,clusterSize);
  }
  //------------------------------
  public void forward(float[] in) {
    //Sy = new float[outputs]; // reset output states
    setAR(Sy,-1);

    clusterManager.forward(in);
    for (int i=0; i<cells.length; i++) {
      if(clusterManager.isActive(i))
        Sy[i] = cells[i].forward(in);
    }
  }
  //------------------------------
  public void backward(float[] dy) {
  //  setAR(dSx,0);
  dSx = new float[dSx.length];
    float[] tmp = new float[dSx.length];
    for (int i=0; i<cells.length; i++) {
      if(clusterManager.isActive(i)){
        tmp = cells[i].backward(dy[i]);
        for (int j=0; j<tmp.length; j++) {
          dSx[j]+=tmp[j];
        }
      }
    }
    // update Cluster Manager:
    clusterManager.backward(dy);
  }
  //------------------------------
  public void update(float lr) {
    //lr*=log(layerLevel+1); // higher learning rate near outputs
    //lr*=1f/layerLevel; // higher learning rate near inputs
    //lr *= layerLevel;
    for (int i=0; i<cells.length; i++)
    cells[i].update(lr);
    clusterManager.update(lr);
  }
  //------------------------------
  public void resetStates() {
    for (int i=0; i<cells.length; i++)
    cells[i].resetStates();
  }
  //------------------------------
  public void setAR(float[] ar,float val){
    for(int i=0;i<ar.length;i++)ar[i]=val;
  }
}

//================================
class MemoryCell {
  float[][] x;
  float Sx;
  float Sh;
  float dSh;
  float dSx;
  float[] dx;
  float[][] Wx, dWx;
  float Bx, dBx;
  int numInputs;
  float weightSize   = 0.10f;
  boolean generalisation = false;
  int steps, last;
  float ret;
  int level;
  //================================
  MemoryCell(int num, int t,int lvl) {
    if (num<1) num=1;
    level = lvl;
    steps     = t;
    last      = t-1;
    numInputs = num;
    weightSize = .1f; // 1.0 / numInputs;
    initArrays();
  }
  //================================
  public void initArrays() {

    x     = new float[steps][numInputs];
    dx    = new float[numInputs];
    Wx    = new float[steps][numInputs];
    dWx   = new float[steps][numInputs];
    for (int t=0; t<steps; t++) {
      for (int i=0; i<numInputs; i++) {
        Wx[t][i] =initWeight(1f/level);
      }
    }
  }
  //================================
  //
  public float forward(float[] _x) {
    push(x, _x);
    //  orig:
    Sx = 0;
    for (int t = 0; t<steps; t++) {
      for (int in=0; in<numInputs; in++) {
        Sx += x[t][in] * Wx[t][in];
      }
    }
    Sx += Bx;
    Sh = softSign(Sx);
    //Sh = Sx;
    return Sh;
  }
  //================================
  public float[] backward(float error) {
    dx = new float[dx.length];

    //dSx = Sx * (1 - Sx ) * error;       // sigmoid
    // dSx = deSoftSign(Sx) * error;         // Softsign
    dSx = (2.0f / sq(1.0f + abs(Sx)))*error;
    //dSx = (1 - sq(Sx)) * error;       // tanh

    //   dSx = Sh==1 ? error:-error;
  //dSx = error;

    dBx += dSx;
    for (int t=steps-1; t>=0; t--) {
      for (int in=0; in<numInputs; in++) {
        dx[in]     += Wx[t][in]*dSx;
        dWx[t][in] += x[t][in] *dSx;
      }
    }
    return dx;
  }
  //================================
  public void update(float learningRate) {
    for (int t=0; t<steps; t++) {
      for (int i=0; i<Wx[t].length; i++) {
        Wx[t][i]  += dWx[t][i]  * learningRate;
      }
    }
    Bx += dBx * learningRate;
    // generalize:
    for (int t=0; t<steps; t++) {
      for (int i=0; i<Wx[t].length; i++) {
        Wx[t][i]  *=(1-learningRate*0.01f);
      }
    }

    resetGradients();
    resetStates(); // why?
  }
  //================================
  private void resetGradients() {
    dWx = new float[dWx.length][dWx[0].length];
    dBx = 0;
  }
  //================================
  public void resetStates() {
    x = new float[steps][numInputs];
  }
  //================================
  // pushes value into array (from top) (bottom-value drops out)
  private void push(float[] ar, float f) {
    for (int i=1; i<ar.length; i++)
    ar[i-1] = ar[i];
    ar[ar.length-1] = f;
  }
  //================================
  private void push(float[][] ar, float[] f) {
    for (int i=1; i<ar.length; i++)
    for (int j=0; j<ar[i].length; j++)
    ar[i-1][j] = ar[i][j];
    for (int j = 0; j < f.length; j++)
    ar[ar.length-1][j] = f[j];
  }
  //================================
  private float sigmoid(float val) {
    return 1 / (1 + exp(-1 * val));
  }
  //================================
  private float tanh(float x) {
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
  }
  //================================
  private float rectifiedLinear(float x) {
    return max(x, 0.0f);
  }
  //================================
  private float softSign(float x) {
    return 2.0f * x / (1.0f + abs(x));
  }
  private float deSoftSign(float x){
    return (2.0f / sq(1.0f + abs(x)));
  }
  //================================
  private float coinflip() {
    if (random(100)>=50) return 1;
    else return -1;
  }
  //================================
  private float initWeight(float w) {
    return random(-w, w);
  }
  //================================
  public float averageWeightSize() {
    return absSum(Wx) / (steps*numInputs);
  }
  //================================
  private float absSum(float[] ar) {
    ret=0;
    for (int i=0; i<ar.length; i++)
    ret+=abs(ar[i]);
    return ret;
  }
  //================================
  private float absSum(float[][] ar) {
    ret = 0;
    for (int i=0; i<ar.length; i++)
    for (int j=0; j<ar[i].length; j++)
    ret += abs(ar[i][j]);
    return ret;
  }
  //================================
}
//------------------------------
class Network {
  Layer[] layers;
  float lr = 1e-2f;
  int batchSize = 0; // if !0 -> batched update every x steps while training
  int cellsPerCluster = 2;
  int memory = 1;
  boolean squaredError = true;
  //------------------------------
  // init
  Network() {
    layers = new Layer[0];
  }
  //------------------------------
  public void resetStates(){
    for(int i=0;i<layers.length;i++)
    layers[i].resetStates();
  }
  //------------------------------
  // learning function
  // recieves dataset with inputs and targets
  // forward -> backprop -> upate
  // returns avg error per output
  public float learn(float[][] in, float[][] out) {
    float[] ret = new float[out[0].length];
    float[] err = new float[out[0].length];

    float absErr = 0;
    float thisErr=0;
    for (int i=0; i<in.length; i++) { // for every entry in the dataset
      ret = forward(in[i]);
      thisErr = 0;
      for (int j=0; j<ret.length; j++) {  // calculate error for this entry
        err[j] = out[i][j]-ret[j];
        thisErr += abs(err[j]);
        if(squaredError) err[j]*=abs(err[j]); // error squared!

      }
      backward(err); // gradient descent
      absErr+=thisErr/ret.length;

      if(batchSize!=0 && i%batchSize==0)
        update(batchSize); // batchlearn
    }

    update(in.length); // update network

    return absErr/in.length;
  }

  //------------------------------
  public float[] forward(float[] in) {
    layers[0].forward(in);
    for (int i=1; i<layers.length; i++) {
      layers[i].forward(layers[i-1].Sy);
    }
    return layers[layers.length-1].Sy;
  }
  //------------------------------
  public void backward(float[] err) {
    layers[layers.length-1].backward(err);
    for (int i=layers.length-2; i>=0; i--) {
      layers[i].backward(layers[i+1].dSx);
    }
  }
  //------------------------------
  public void update(int numSamples) {
    for (int i=0; i<layers.length; i++)
    layers[i].update(lr/numSamples);
  }
  //------------------------------
  public void addLayer(int _in, int _out) {
    layers =  (Layer[]) expand(layers, layers.length+1);
    layers[layers.length-1] = new Layer(_in, _out, cellsPerCluster,memory,layers.length);
    println("adding layer: "+_in+" -> "+_out);
  }
  //------------------------------
}
class ClusterManager {

  float[]w,dw;
  float[] b,db;
  float[] Sh;
  int clusterSize;
  boolean[] Sy;
  int numClusters;
  //float[] S0;
  float[] S0;
  int connectedNeurons;
  //------------------------------
  // setup
  ClusterManager(int in, int out, int clusterSize) {
    this.clusterSize=clusterSize;
    this.connectedNeurons = out;
    numClusters = (out / clusterSize>0?out/clusterSize:1 );
    numClusters+=clusterSize >= out ? 0:1;

    S0 = new float[in];
    Sy = new boolean[numClusters];
    w = new float[numClusters*in];
    dw = new float[numClusters*in];
    b = new float[numClusters];
    db = new float[numClusters];
    Sh = new float[numClusters];
    for(int i=0;i<w.length;i++){
      w[i]=random(-1,1)*0.02f;
    }
    //println("c: "+numClusters);

  }
  //------------------------------
  public float sum(float ar[]){
    float ret = 0;
    for(int i=0;i<ar.length;i++)
    ret+=(ar[i]);
    return ret;
  }
  //------------------------------
  public void forward(float[] in) {

    S0 = in;
    Sh = new float[Sh.length];


    //    Sy[0] = true; // force one cluster open at all times
    //-> remove backprop for this aswell...^^
    // shouldn't really be necessary, but could speed up the early learning phase..

    for (int i=0; i<numClusters; i++) {
      for(int j=0;j<in.length;j++){
        if(i*clusterSize+j>=w.length)continue;
        Sh[i] += in[j] * w[i*clusterSize+j];
      }
      Sh[i]+=b[i];
      // add softsign!?
      //Sh[i] = softSign(Sh[i]);
      Sy[i] = Sh[i]>.1f ? true:false;
    }

  }
  //------------------------------
  public void backward(float[] dy) {

    // dy is of length clusterSize*numClusters! -> of length outputs

    for (int i=0; i<numClusters; i++) {

      for(int k=0;k<clusterSize;k++){
        // with this the weights stay near the on/off state...// makes sense? idk...
        //  if(Sy[i]&&  dy[i*clusterSize+k]>0 || !Sy[i] &&  dy[i*clusterSize+k]<0) continue;
          if(i*clusterSize+k>=dy.length)continue;
        //dy[i*clusterSize+k] = deSoftSign(Sh[i])*dy[i*clusterSize+k];
        db[i]+= dy[i*clusterSize+k];

        for(int j=0;j<S0.length;j++){
          if(i*clusterSize+j>=w.length)continue;
          dw[i*clusterSize+j] =+ S0[j] * dy[i*clusterSize+k];
        }
      }
    }
  }
  //------------------------------
  public void update(float lr) {
    //  lr*=1e-1; // clusterManager learns slower than actual NN
    for (int i=0; i<b.length; i++){
      b[i]+=db[i]*lr;
      db[i]=0;
    }
    for(int i=0;i<w.length;i++){
      w[i]+=dw[i]*lr;
      //w[i]*=(1f-lr*1e-2); // gerneralisation -> does this make any sense? idk...
      dw[i]=0;
    }
  }
  //------------------------------
  public boolean isActive(int index){
    return Sy[inputToCluster(index)];   // sth like this!
  }
  //------------------------------
  public int inputToCluster(int index){
    return index/clusterSize < numClusters ?
      index/clusterSize:
      numClusters;
  }
  //------------------------------

  private float softSign(float x) {
    return 2.0f * x / (2.0f + abs(x));
  }
  //------------------------------
  private float deSoftSign(float x){
    return (2.0f / sq(2.0f + abs(x)));
  }
  //------------------------------
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "--full-screen", "--bgcolor=#050505", "--hide-stop", "sketch_180727a_clusteredNetwork" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
