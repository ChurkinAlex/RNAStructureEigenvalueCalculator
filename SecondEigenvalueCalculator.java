import java.io.BufferedReader;


import Jama.Matrix;
import Jama.EigenvalueDecomposition;
import java.io.FileReader;
import java.io.IOException;

public class SecondEigenvalueCalculator {
	
	public static double calcEig(double[][] mat) {
	      Matrix A = new Matrix(mat);
	      EigenvalueDecomposition e = A.eig();
	      double [] realPart = e.getRealEigenvalues();
	      
	      if (realPart.length<2) return 0;
	      return realPart[1];
	}
	
	public static double[][] makeLaplacianMatrix(String shapiro){
		  double[][] matrix;
		  int sizeOfMatrix , sizeOfShapiro;
		  int numOfParens=0;
		  int node1=0, node2=0;
		  int isNode=0;
		  int isS=0;
		  int count;
		  String newShapiro="";
		  
		  sizeOfShapiro = shapiro.length();
		  sizeOfMatrix = 0;
		   
		  for (int i=0; i<sizeOfShapiro; i++){
		    if ((((int) shapiro.charAt(i))>64 && ((int) shapiro.charAt(i))<91
		         && shapiro.charAt(i)!='S' && shapiro.charAt(i)!='E')&&
		        !((shapiro.charAt(i)=='I') && (shapiro.charAt(i+1)=='2') && (shapiro.charAt(i+2)==')'))&&
		        !((shapiro.charAt(i)=='B') && (shapiro.charAt(i+1)=='1') && (shapiro.charAt(i+2)==')')))
		      sizeOfMatrix++;
		  }
		  
		  matrix = new double[sizeOfMatrix][sizeOfMatrix];
		  
		  count = -1;
		  for (int i=0; i<sizeOfShapiro; i++){
		    if (shapiro.charAt(i) == 'S') isS=1;
		    else if ((shapiro.charAt(i)=='I') && (shapiro.charAt(i+1)=='2') && (shapiro.charAt(i+2)==')')){
		      newShapiro+= 'I';
		    }
		    else if ((shapiro.charAt(i)=='B') && (shapiro.charAt(i+1)=='1') && (shapiro.charAt(i+2)==')')){
		      newShapiro+= 'B';
		    }
		    else if (((int) shapiro.charAt(i))>64 && ((int) shapiro.charAt(i))<91){
		      count++;
		      isNode = 1;
		    }
		    else if (((int) shapiro.charAt(i))>47 && ((int) shapiro.charAt(i))<58){
		    }
		    else if (isNode!=0){
		      isNode = 0;
		      newShapiro+= (char) ((count / 10)+48);
		      newShapiro+= (char) ((count % 10)+48);
		      i--;
		    }
		    else if (isS!=0){
		      isS = 0;
		      newShapiro+= 'S';
		      i--;
		    }
		    else{
		      newShapiro+= shapiro.charAt(i);
		    }
		  }

		  //System.out.println(newShapiro);
		  for (int i=0; i<newShapiro.length(); i++){
		    if (((int) newShapiro.charAt(i))>47 && ((int) newShapiro.charAt(i))<58){
		      node1 = 10*(((int) newShapiro.charAt(i)) - 48) + (((int) newShapiro.charAt(i+1)) - 48);
		      i = i+5;
		      numOfParens = 0;
		      for (int k=i; k<newShapiro.length(); k++){
		        if (newShapiro.charAt(k) == '(') numOfParens ++;
		        if (newShapiro.charAt(k) == ')') numOfParens --;
		        if ((newShapiro.charAt(k) == 'I' || newShapiro.charAt(k) == 'B') && numOfParens <=0) {
		          numOfParens=0;
		          k = k+3;
		        }
		        if (((int) newShapiro.charAt(k))>47 &&
		            ((int) newShapiro.charAt(k))<58 && numOfParens==0){
		        node2 =  10*(((int) newShapiro.charAt(k)) - 48) + (((int) newShapiro.charAt(k+1)) - 48);
		          matrix[node1][node1]++;
		          matrix[node2][node2]++;
		          matrix[node1][node2]=-1;
		          matrix[node2][node1]=-1;
		          break;
		        }
		      }
		      i--;
		    }
		  }
		  
		  return matrix;
		}

	
	public static String reverseComplement(String seq) {
		String rev = "";
		for(int i=seq.length()-1; i>=0; i--) {
			char c = seq.charAt(i);
			char rc = '*';
			if(c=='A' || c=='a') rc = 'T';
			else if (c=='T' || c=='t' || c=='U' || c=='u') rc = 'A';
			else if (c=='C' || c=='c') rc = 'G';
			else if (c=='G' || c=='g') rc = 'C';
			else if (c=='Y' || c=='y') rc = 'R';
			else if (c=='R' || c=='r') rc = 'Y';
			else if (c=='S' || c=='s') rc = 'S';
			else if (c=='W' || c=='w') rc = 'W';
			else if (c=='K' || c=='k') rc = 'M';
			else if (c=='M' || c=='m') rc = 'K';
			else if (c=='B' || c=='b') rc = 'V';
			else if (c=='V' || c=='v') rc = 'B';
			else if (c=='D' || c=='d') rc = 'H';
			else if (c=='H' || c=='h') rc = 'D';
			else if (c=='N' || c=='n') rc = 'N';
			rev+=rc;
		}
		
		return rev;
	}
	
	public static void analyzeStructs() throws IOException{
		String file = "D:\\structs2.txt";

		BufferedReader bufferedReader = new BufferedReader(new FileReader(file));
        String curLine;

        String seq = "";
        int count = 0;
        double max = 0;
        double min = 10;
        int linear = 0;
        int non_linear = 0;
        while ((curLine = bufferedReader.readLine()) != null){
        	//if (curLine.contains("Undifined")) {
        	if (curLine.contains(">")) {
        		seq = bufferedReader.readLine();
        		String dot = bufferedReader.readLine();
        		String shap = bufferedReader.readLine();
        		bufferedReader.readLine();
        		double eig = calcEig(makeLaplacianMatrix(shap));
        		if(eig>max)
        			max = eig;
        		if(eig<min)
        			min = eig;
        		//GAAC

        		int M=0, H=0;
        		for(int i=0; i<shap.length(); i++) {
        			if (shap.charAt(i)=='H') H++;
        			if (shap.charAt(i)=='M') M++;
        		}
        		/*if (H>2) {
        			if(seq.contains("GAAC")) {
        				for(int i=0; i<seq.length()-3; i++) {
        					if(seq.charAt(i)=='G' && seq.charAt(i+1)=='A' && seq.charAt(i+2)=='A' && seq.charAt(i+3)=='C') {
        						if (dot.charAt(i)=='.' && dot.charAt(i+1)=='.' && dot.charAt(i+2)=='.' && dot.charAt(i+3)=='.') {
        							if (i+6 >= seq.length()-3) continue;
        							System.out.println(curLine);
        							System.out.println(seq.subSequence(i-3, i+6));
        							System.out.println(dot.subSequence(i-3, i+6));
        						}
        					}
        				}
        			}
        		}*/
        		
        		/*System.out.println(curLine);
        		System.out.println(shap+"\n"+eig);
        		System.out.println("H="+H+" M="+M+"\n");*/
        		if (H==2) linear++;
        		else non_linear++;
        	}
        }
        System.out.println(min+"\n"+max+"\n"+linear+"\n"+non_linear);
        bufferedReader.close();
	}
	
	
	
	public static void main(String[] args) throws IOException  {
		analyzeStructs();		
	}

	
	
	
	
	
	
	
    
    
    


}
