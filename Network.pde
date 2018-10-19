//------------------------------
class Network {
  Layer[] layers;
  float lr = 1e-2;
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
  void resetStates(){
    for(int i=0;i<layers.length;i++)
    layers[i].resetStates();
  }
  //------------------------------
  // learning function
  // recieves dataset with inputs and targets
  // forward -> backprop -> upate
  // returns avg error per output
  float learn(float[][] in, float[][] out) {
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
  float[] forward(float[] in) {
    layers[0].forward(in);
    for (int i=1; i<layers.length; i++) {
      layers[i].forward(layers[i-1].Sy);
    }
    return layers[layers.length-1].Sy;
  }
  //------------------------------
  void backward(float[] err) {
    layers[layers.length-1].backward(err);
    for (int i=layers.length-2; i>=0; i--) {
      layers[i].backward(layers[i+1].dSx);
    }
  }
  //------------------------------
  void update(int numSamples) {
    for (int i=0; i<layers.length; i++)
    layers[i].update(lr/numSamples);
  }
  //------------------------------
  void addLayer(int _in, int _out) {
    layers =  (Layer[]) expand(layers, layers.length+1);
    layers[layers.length-1] = new Layer(_in, _out, cellsPerCluster,memory,layers.length);
    println("adding layer: "+_in+" -> "+_out);
  }
  //------------------------------
}
