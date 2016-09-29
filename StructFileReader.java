import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class StructFileReader {

	HashMap<String, double[]> structureSequences;
	HashMap<String, boolean[]> isValidScore; 
	Double mean;
	Double variance;
	
	public StructFileReader(){
		structureSequences = new HashMap<String, double[]>(); // scores
		isValidScore = new HashMap<String, boolean[]>();	  // keep track of NA scores
	}
	
	public void readFile(String infile, Integer minLength, Double minScore) throws IOException{
		// read the file into two hashmaps and calculate the mean score

		BufferedReader br = new BufferedReader(new FileReader(infile));
		String structLine;
		double sum = 0;
		int total = 0;
		
		while ((structLine = br.readLine()) != null){
			String[] elts = structLine.split(",");
			if (minLength != null && elts.length >= minLength){	
	
				double[] scores = new double[elts.length - 1];
				boolean[] validity = new boolean[elts.length - 1];
				for (int i=0;i<scores.length;i++){
					if (elts[i+1].compareTo("NA") == 0 || elts[i+1].compareTo("NaN") == 0){
						scores[i] = 0;
						validity[i] = false;
					}	
					else {
						double score = Double.parseDouble(elts[i+1]);

						if (minScore != null && Math.abs(score) < minScore){
							scores[i] = 0;
							validity[i] = false;
						}
						
						else {
							scores[i] = score;
							validity[i] = true;
							sum += score;
							total += 1;
						}
					}
				}
				structureSequences.put(elts[0], scores);
				isValidScore.put(elts[0], validity);
			}
		}
		br.close();
		mean = sum/total;
	}

	
	public HashMap<String, double[]> getStructureSequence(){
		return structureSequences;
	}
	
	public HashMap<String, boolean[]> getIsValidScore(){
		return isValidScore;
	}
	
	public Double getMean(){
		return mean;
	}
	
	public Double getVariance(){

		// if the variance is already computed return it
		if (variance != null) return variance;
		
		// otherwise pass through the data and compute it
		double sumOfSquares = 0;
		int total = 0;
		for (String k : structureSequences.keySet()){
			double [] seq = structureSequences.get(k);
			boolean[] val = isValidScore.get(k);
			for (int i=0;i<seq.length;i++){
				if (val[i]){
					sumOfSquares += Math.pow((seq[i] - mean),2);
					total += 1;
				}
			}
		}
		variance = Math.pow((sumOfSquares/total),0.5);
		return variance;
	}
}
