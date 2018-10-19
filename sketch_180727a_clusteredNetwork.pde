import java.util.*;

Network nn = new Network();
float[] testErrorList=new float[0];
float[] activeCells=new float[0];
float minCells = 100;
float err = 0;
int timer = 0;
PFont font;
//------------------------------
void setup() {
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

  nn.lr = 1.0e-1; // learningRate


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

void draw() {
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
          output += nn.layers[l].clusterManager.Sy[j]? "■":"□";
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
void drawLine(float y, float x) {
  strokeWeight(3);
  line(0, y, 0, map(x, 2, -2, -height/3, height/3));
  translate(4, 0);
}
//------------------------------
int highestIndex(float ar[]){
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
void displayError(float[] list) {
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
  void keyPressed() {
    if (key=='+') nn.lr*=2;
    if (key=='-') nn.lr/=2;
  }
