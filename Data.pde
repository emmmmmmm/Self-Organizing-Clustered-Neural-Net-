float[][] inputData, targetData;
float[][] trainingIn,trainingOut,validIn,validOut;

void buildData() {
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
    inputData[i][0] = sin(i*.1);
    // inputData[i][1] = cos(i*.1);
    //}
  }
  for (int i=0; i<inputData.length; i++) {
    targetData[i][0] =cos(i*.2);
  }

  for (int i=1; i<inputData.length; i++) {
    for (int j=0; j<inputData[i].length; j++) {
      inputData[i][0] =sin(i*.1);
      // inputData[i][1] = cos(i*.1);
    }
  }
  for (int i=0; i<inputData.length; i++) {

    targetData[i][0] = abs(sin(i*.1)*cos(i*.1));
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
void importData() {
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
    f[i][0] = int(tmp[i][0]);
    f[i][1] = int(tmp[i][1]);
    f[i][2] = int(tmp[i][2]);
    f[i][3] = int(tmp[i][3]);
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
    f[i][0] = int(tmp[i][0]);
    f[i][1] = int(tmp[i][1]);
    f[i][2] = int(tmp[i][2]);
    f[i][3] = int(tmp[i][3]);
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

void shuffleArray(float[][] ar1, float[][] ar2)
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


void shuffleArray(String[] array)
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
