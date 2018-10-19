/*

//------------------------------
class Layer {
  CellCluster cellCluster[];
  float[] Sy;
  float[] dSx;
  int inputs,  outputs;
  //------------------------------
  Layer(int _in, int _out, int cellClusterPerCluster) {
    inputs=_in;
    outputs = _out;
    Sy = new float[outputs];
    dSx = new float[inputs];
    cellCluster = new CellCluster[outputs];

    for (int i=0; i<outputs; i++) {
      cellCluster[i] = new CellCluster(inputs, cellClusterPerCluster);
    }
  }
  //------------------------------
  void forward(float[] in) {
    for (int i=0; i<cellCluster.length; i++) {
      Sy[i] = cellCluster[i].forward(in);
    }
  }
  //------------------------------
  void backward(float[] dy) {
    dSx = new float[dSx.length];
    float[] tmp = new float[dSx.length];
    for (int i=0; i<cellCluster.length; i++) {
      tmp = cellCluster[i].backward(dy[i]);
      for (int j=0; j<tmp.length; j++) {
        dSx[j]+=tmp[j];
      }
    }
  }
  //------------------------------
  void update(float lr) {
    for (int i=0; i<cellCluster.length; i++)
      cellCluster[i].update(lr);
  }
  //------------------------------
  void resetStates() {
    for (int i=0; i<cellCluster.length; i++)
      cellCluster[i].resetStates();
  }
}
*/
