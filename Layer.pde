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
  void forward(float[] in) {
    //Sy = new float[outputs]; // reset output states
    setAR(Sy,-1);

    clusterManager.forward(in);
    for (int i=0; i<cells.length; i++) {
      if(clusterManager.isActive(i))
        Sy[i] = cells[i].forward(in);
    }
  }
  //------------------------------
  void backward(float[] dy) {
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
  void update(float lr) {
    //lr*=log(layerLevel+1); // higher learning rate near outputs
    //lr*=1f/layerLevel; // higher learning rate near inputs
    //lr *= layerLevel;
    for (int i=0; i<cells.length; i++)
    cells[i].update(lr);
    clusterManager.update(lr);
  }
  //------------------------------
  void resetStates() {
    for (int i=0; i<cells.length; i++)
    cells[i].resetStates();
  }
  //------------------------------
  void setAR(float[] ar,float val){
    for(int i=0;i<ar.length;i++)ar[i]=val;
  }
}
