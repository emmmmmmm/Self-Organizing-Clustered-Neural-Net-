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
      w[i]=random(-1,1*1.2f);
    }
    //println("c: "+numClusters);

  }
  //------------------------------
  float sum(float ar[]){
    float ret = 0;
    for(int i=0;i<ar.length;i++)
    ret+=(ar[i]);
    return ret;
  }
  //------------------------------
  void forward(float[] in) {

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
      Sy[i] = Sh[i]>.1 ? true:false;
    }

  }
  //------------------------------
  void backward(float[] dy) {

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
  void update(float lr) {
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
  boolean isActive(int index){
    return Sy[inputToCluster(index)];   // sth like this!
  }
  //------------------------------
  int inputToCluster(int index){
    return index/clusterSize < numClusters ?
      index/clusterSize:
      numClusters;
  }
  //------------------------------

  private float softSign(float x) {
    return 2.0 * x / (2.0 + abs(x));
  }
  //------------------------------
  private float deSoftSign(float x){
    return (2.0 / sq(2.0 + abs(x)));
  }
  //------------------------------
}
