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
