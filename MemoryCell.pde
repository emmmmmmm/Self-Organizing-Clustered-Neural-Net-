
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
  float weightSize   = 0.10;
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
    weightSize = .1; // 1.0 / numInputs;
    initArrays();
  }
  //================================
  void initArrays() {

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
    dSx = (2.0 / sq(1.0 + abs(Sx)))*error;
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
        Wx[t][i]  *=(1-learningRate*0.01);
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
    return max(x, 0.0);
  }
  //================================
  private float softSign(float x) {
    return 2.0 * x / (1.0 + abs(x));
  }
  private float deSoftSign(float x){
    return (2.0 / sq(1.0 + abs(x)));
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
