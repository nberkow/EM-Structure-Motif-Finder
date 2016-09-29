import java.util.HashMap;
import java.util.Random;
import java.util.Set;

public class EMMotifFinder {

	HashMap<String, Integer> startPositions;
	HashMap<String, Double> startPositionProbabilities;
	double backgroundDistMean;
	double backgroundDistVariance;
	double [] mu;
	double [] sigma;
	
	// The difference between iterations
	double delta;
	
	// The difference between iterations where they are considered to have converged
	final double EPSILON = Math.pow(10, -200);
	
	// The number of positions in the motif and how spaced out they are
	int motifPositions;
	int motifSpacing;
	int motifLength;
	
	HashMap<String, double []> sequences; // Track the sequence of scores
	HashMap<String, boolean []> validity; // Track where NA values occurred in the file
	Random rgen = new Random();
	int totalIterations = 0;
	Set<String> seqNames;
	
	// These are for visualizing the motif after the algorithm is run
	HashMap<String, double[]> motifCenteredScoreWindow;
	HashMap<String, boolean []> motifCenteredValidityWindow;

	public EMMotifFinder (int mPos, int mSpace, HashMap<String, double []> seq, HashMap<String, boolean []> val){
		motifPositions = mPos;
		motifSpacing = mSpace;
		motifLength = motifPositions * (motifSpacing) + 1;
		
		sequences = seq;
		validity = val;
		seqNames = seq.keySet();
		startPositions = new HashMap<String, Integer>();
		startPositionProbabilities = new HashMap<String, Double>();
		for (String k : seqNames){
			startPositions.put(k, -1);
		}
	}
	
	// Set the initial values of mu and sigma
	public void initializeStartingPoints(double [] m, double [] s, boolean randomize, double bgm, double bgv){
		delta = Math.pow(10, 200);
		mu = m;
		sigma = s;
		backgroundDistMean = bgm;
		backgroundDistVariance = bgv;
	}

	// Set the initial values of mu and sigma
	public void initializeStartingPoints(double m, double s, boolean randomize){
		delta = Math.pow(10, 200);
		mu = new double[motifPositions];
		sigma = new double[motifPositions];

		for (int i=0; i<motifPositions;i++){
			if (randomize){
				mu[i] = (rgen.nextGaussian() * s) + m;
			}
			else {
				mu[i] = m;
			}
			sigma[i] = s;
		}
	}

	// This version of the eStep maximizes a probability value

	/*//E Step - maximize the motif starting positions using mu and sigma
	private void eStep(){
		// Run the EM motif find on a real dataset. 
		// Calculate probabilities for all possible motif starting positions
		// 				and
		// Find the position with maximum probability for each sequence and update delta

		delta = 0;
		for (String k : seqNames){

			// record the start position from the previous iteration
			int previousStart;
			if (startPositions.containsKey(k)) previousStart = startPositions.get(k);
			else previousStart = -100;
			startPositions.put(k, 0);

			// check the probability of each position and test whether it is the max
			double[]seq = sequences.get(k);
			boolean[]val = validity.get(k);
			double maxP = -1;
			
			for (int i=0;i<=seq.length-(motifLength)+1;i++){
				double p = 1;
				for (int j=0;j<motifPositions;j++){
					
					// If the motif contains an NA, the motif has a probability of 0
					if (!val[i + (j * motifSpacing)]) p = 0; 
					p = p * dnorm(seq[i + (j * motifSpacing)], mu[j], sigma[j]) + EPSILON;
					//System.out.println(">>\t" + i + "\t" + j + "\t" + 
					//seq[i + (j * motifSpacing)] + "\t" + mu[j] + "\t" 
					//+ sigma[j] + "\t" + maxP + "\t" + p + "\t" + dnorm(seq[i + (j * motifSpacing)], mu[j], sigma[j]) );
				}

				if (p > maxP){
					startPositions.put(k,  i);
					startPositionProbabilities.put(k, p);
					maxP = p;
				}
			}
			
			//System.out.println(">\t" + k + "\t" + maxP + "\t" + startPositions.get(k));
			
			delta += Math.abs(previousStart - startPositions.get(k));
		}	
	}*/
	

	//E Step - maximize the motif starting positions using mu and sigma
	private void eStep(){
		// Run the EM motif find on a real dataset. 
		// Calculate probabilities for all possible motif starting positions
		// 				and
		// Find the position with maximum probability for each sequence and update delta

		delta = 0;
		for (String k : seqNames){

			// record the start position from the previous iteration
			int previousStart;
			if (startPositions.containsKey(k)) previousStart = startPositions.get(k);
			else previousStart = -100;
			startPositions.put(k, 0);

			// check the probability of each position and test whether it is the max
			double[]seq = sequences.get(k);
			boolean[]val = validity.get(k);
			double min_distance = 100;
			
			for (int i=0;i<=seq.length-(motifLength)+1;i++){
				double dist = 0;
				for (int j=0;j<motifPositions;j++){
					
					// If the motif contains an NA, the motif has a probability of 0
					if (!val[i + (j * motifSpacing)]) p = dist = 100;
 
					// Othewise add up the geometric distances of all the scores
					dist = dist + (seq[i + (j * motifSpacing)] - mu[j])**2
					dist = dist**0.5

					//p = p * xdnorm(seq[i + (j * motifSpacing)], mu[j], sigma[j]) + EPSILON;
					//System.out.println(">>\t" + i + "\t" + j + "\t" + 
					//seq[i + (j * motifSpacing)] + "\t" + mu[j] + "\t" 
					//+ sigma[j] + "\t" + maxP + "\t" + p + "\t" + dnorm(seq[i + (j * motifSpacing)], mu[j], sigma[j]) );
				}

				if (dist < min_distance){
					startPositions.put(k,  i);
					startPositionProbabilities.put(k, p);
					min_distance = dist;
				}
			}
			
			//System.out.println(">\t" + k + "\t" + maxP + "\t" + startPositions.get(k));
			
			delta += Math.abs(previousStart - startPositions.get(k));
		}	
	}
	


	/*
	private void debug(String seqName, int pos, int i, int j, double[] seq, int startPos, double [][] scoreMatrix){
		System.out.println(seqName);
		System.out.println("i:\t" + i + "\n");
		System.out.println("j:\t" + j + "\n");
		System.out.println("pos:\t" + pos + "\n");
		System.out.println("startPos:\t" + startPos + "\n");
		System.out.println("seqLen:\t" + seq.length + "\n");
		System.out.println("matrixLen:\t" + scoreMatrix.length + "\t" + scoreMatrix[0].length + "\n");
	}*/

	//M Step - compute mu and sigma using the motif starting positions
	private void mStep(){
		
		// Find the weight contribution of each position contributing to the motif
		// Keep track of the row and column sums for scaling
		double [][] weightMatrix = new double [sequences.size()][motifPositions];
		double [][] scoreMatrix = new double [sequences.size()][motifPositions];
		double [] rowSums = new double [sequences.size()];
		double [] colSums = new double [motifPositions];
		int i = 0;
		for (String k : seqNames){
			double [] seq = sequences.get(k);
			boolean [] val = validity.get(k);
			
			int startPos = startPositions.get(k);
			for (int j=0; j<motifPositions;j++){	
				int pos = startPos + (j * motifSpacing);
				
				double score = seq[pos];
				boolean v = val[pos];
				scoreMatrix[i][j] = score;
				
				double wc;
				if (v)	 wc = dnorm(score, mu[j], sigma[j]);
				else wc = EPSILON;
				
				weightMatrix[i][j] = wc;
				rowSums[i] += wc;
				colSums[j] += wc;
			}
			i ++;
		}		
		
		// Calculate the weight of each sequence
		double [] sequenceWeights = new double [sequences.size()];
		double weightSum = EPSILON;
		for (int x=0;x<sequences.size();x++){
			sequenceWeights[x] = 1;
			for (int y=0;y<motifPositions;y++){
				sequenceWeights[x] = sequenceWeights[x] * weightMatrix[x][y] / (colSums[y]) + EPSILON;
			}
			weightSum += sequenceWeights[x];
		}
		
		// re-scale the weights so that their mean is one
		for (int w=0;w<sequenceWeights.length;w++){
			sequenceWeights[w] = sequences.size() * sequenceWeights[w] / (weightSum);	
		}

		// calculate the weighted means
		for (int y=0;y<motifPositions;y++){
			mu[y] = 0;
			for (int x=0;x<sequences.size();x++){
				mu[y] += scoreMatrix[x][y] * sequenceWeights[x] / sequences.size();
			}
		}

		// calculate sum of squares using un-weighted scores
		double [] squaredDiffFromMean = new double[motifPositions];
		for (int y=0; y<motifPositions;y++){
			squaredDiffFromMean[y] = 0;
			for (int x=0;x<sequences.size();x++){
				squaredDiffFromMean[y] += Math.pow(mu[y] - scoreMatrix[x][y], 2);
			}
		}

		// calculate sigma
		for (int x=0;x<motifPositions;x++){
			sigma[x] = Math.pow(squaredDiffFromMean[x]/sequences.size(), 0.5) + .001;
			if (sigma[x] < backgroundDistVariance/2){
				sigma[x] = backgroundDistVariance/2;
			}
		}
	}

	public void runEM(int maxIterations, double deltaMin){
		
		// run the algorithm to completion
		while(totalIterations <= maxIterations && delta > deltaMin){

			// Do one iteration of EM	
			eStep();
			mStep();
			totalIterations += 1;		
		}
	}

	// returns the height of the PDF for x. Identical to dnorm function in R
	public double dnorm(double x, double m, double s){
		double p = 1 - EPSILON;
		if (s != 0){
			double q1 = 1/Math.pow((2 * Math.PI * Math.pow(s,2)),0.5);
			double q2 = -1 * (Math.pow((x - m),2))/((2 * Math.pow(s,2)));
			p = q1 * Math.pow(Math.E, q2);
		}
		return p + EPSILON;
	}

	public String getPrintableParameters(){
		String printable = "mu:\t";
		for (double m : mu){
			printable = printable + "\t" + m;
		}
		printable += "\nsigma:\t";

		for (double s : sigma){
			printable = printable + "\t" + s;
		}

		printable = printable + "\ndelta:\t" + delta; 

		return printable;
	}

	public String getPrintableStartingPoints(){
		String printable = "";
		for (String k : startPositions.keySet()){
			printable = printable + k + "\t" + startPositions.get(k) + "\t" + startPositionProbabilities.get(k) + "\n";
		}
		return printable;
	}
	public HashMap<String, Integer> getStartPositions() {
		return startPositions;
	}

	public HashMap<String, Double> getStartPosProbabilities() {
		return startPositionProbabilities;
	}
	
	public Double getStartPosSetLikelihood() {
		double p = 1.0;
		for (String s : startPositionProbabilities.keySet()){
			p = p * (startPositionProbabilities.get(s) + EPSILON);
		}
		return p;
	}

	public Set<String> getSeqNames(){
		return seqNames;
	}
	public double getLikelihoodOfMu() {
		double p = 1.0;
		for (double m : mu){
			p = p * dnorm(m, backgroundDistMean, backgroundDistVariance);
		}
		return p;
	}
	public Integer getMotifLength() {
		return motifLength;
	}
	
	public void findMotifCenteredWinows(int flank){
		motifCenteredScoreWindow = new HashMap<String, double[]>();
		motifCenteredValidityWindow = new HashMap<String, boolean[]>();
		
		for (String seqName : startPositions.keySet()){
			double[] fullProfile = sequences.get(seqName);
			boolean[] fullValidity = validity.get(seqName);
			int startPos = startPositions.get(seqName);
			
			int windowStart = startPos - flank;
			int windowEnd = startPos + motifLength + flank;
			int windowSize = windowEnd - windowStart - 1;
			double[] scoreWindow = new double[windowSize];
			boolean[] valWindow = new boolean[windowSize];
			
			for (int i = 0; i < windowSize; i++){
				int pos = windowStart + i;
				if (pos < 0 || pos >= fullProfile.length){
					scoreWindow[i] = 0.0;
					valWindow[i] = false;
				}
				else {
					scoreWindow[i] = fullProfile[pos];
					valWindow[i] = fullValidity[pos];
				}
			}
			motifCenteredScoreWindow.put(seqName, scoreWindow);
			motifCenteredValidityWindow.put(seqName, valWindow);
		}
	}
	
	public HashMap<String, double[]> getScoreWindow(){
		return motifCenteredScoreWindow;
	}
	
	public HashMap<String, boolean[]> getValidityWindow(){
		return motifCenteredValidityWindow;
	}
}
