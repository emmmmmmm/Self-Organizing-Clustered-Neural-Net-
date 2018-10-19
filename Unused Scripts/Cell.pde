class SimpleCell {
  float[] Wx;
  float Bx;
  float Sy;
  float[] Sx;
  float[] dWx;
  float dBx;
  float[] dSx;
  int inputs;
  int i;

  //------------------------------
  SimpleCell(int in) {
    inputs = in;

    Wx= new float[in];
    dWx = new float[in];
    Sx = new float[in];
    dSx = new float[in];
    for (i=0; i<Wx.length; i++)
      Wx[i] = random(-1,1);
  }
  //------------------------------
  float forward(float[] in) {

    Sx = in;
    Sy = 0;
    for (i=0; i<in.length; i++)
      Sy+=Wx[i]*in[i];
    Sy+=Bx;
    Sy= softSign(Sy);
    return Sy;
  }
  //------------------------------
  float[] backward(float dSy) {

    dSy = deSoftSign(Sy) * dSy; // softsign
    dBx += dSy;
    setAR(dSx,0);
    for (i=0; i<Wx.length; i++) {
      dWx[i] += Sx[i]*dSy;
      dSx[i] += Wx[i]*dSy;
    }
    return dSx;
  }
  //------------------------------
  void update(float lr) {
    for (i=0; i<Wx.length; i++) {
      Wx[i]+=dWx[i]*(lr);
    }

    Bx+=dBx*lr;
    setAR(dWx,0);
    dBx=0;
  }
  void resetStates() {
  }
  //------------------------------
  float softSign(float x) {
    return 2.0 * x / (1.0 + abs(x));
  }
  //------------------------------
  float deSoftSign(float y) {
    return (2.0 / sq(1.0+abs(y)));
  }
   //------------------------------
    void setAR(float[] ar,float val){
      for(int i=0;i<ar.length;i++)ar[i]=val;
    }
}
