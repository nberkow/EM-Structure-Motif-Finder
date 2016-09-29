import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;


/**
 * @author nberk
 *
 * A creates and runs EMMotifFinder objects 
 */
public class EM_Manager {
	
	final double EPSILON = Double.MIN_VALUE;
	double min_motif_diff = 20;
	
	List <EMMotifFinder> models = new ArrayList<EMMotifFinder>();
	
	public void runWithRandomStarts(	int iterations, int motifLength, int motifSpacing, 
										double mu, double sigma, 
										HashMap<String, double []> sequences, HashMap<String, boolean []>validity){
		
		for (int i=0;i<iterations;i++){			
			EMMotifFinder emm = new EMMotifFinder(motifLength, motifSpacing, sequences, validity);
			emm.initializeStartingPoints(mu, sigma, true);
			emm.runEM(100, .0001);
			emm.findMotifCenteredWinows(150);
			models.add(emm);
			validity = eraseMotifs(validity, emm.getStartPositions(), emm.getMotifLength());
		}
	}

	/* 	Run the motif finder from specified starting values of mu. Each element
	 	of mu[] is used to make the starting mu for a single run. Currently this 
	 	supports flat starting mu's.	 	 L
	*/
	public void runWithFixedStarts(		int motifLength, int motifSpacing, 
										double [] mu, double [] sigma, 
										HashMap<String, double []>sequences, HashMap<String, boolean []>validity,
										boolean randomize){
		
		for (int i=0;i<mu.length;i++){			
			EMMotifFinder emm = new EMMotifFinder(motifLength, motifSpacing, sequences, validity);
			emm.initializeStartingPoints(mu[i], sigma[i], randomize);
			emm.runEM(100, .0001);
			models.add(emm);
		}
	}
	
	/* 	Run the motif finder from optimized starting mu's. The optimization step
	 *  runs one iteration of the EM Motif Finder for every possible substring
	 *  in the input sequences. These are then used to recompute mu. The mu values
	 *  that are most different from the background distribution are then used
	 *  to run the algorithm to completion.
	 */
	public void runWithOptimizedStarts(	int motifPositions, int motifSpacing, int numberOfStarts, double mmd,
										int maxIterations, double deltaMin,
										double sigmaVal, double bgMean, double bgVar, 
										HashMap<String, double []> sequences, 
										HashMap<String, boolean []>validity){
		
		// get a list of ranked, unique starting motifs
		min_motif_diff = mmd;
		List<double[]> subsequences = getUniqueSubsequences(sequences, motifPositions, motifSpacing, motifPositions * motifSpacing);
		List<double[]> rankedSubsequences = getRankedSubsequences(subsequences, motifPositions, motifSpacing, sigmaVal, bgMean, bgVar, sequences, validity);

		// intialize the starting sigmas
		double [] sigmas = new double[motifPositions];
		for (int s=0;s<motifPositions;s++) sigmas[s] = sigmaVal;
		
		for (int m=0;m<numberOfStarts;m++){
			System.out.println("running EM from starting motif " + m);
			System.out.println(sequenceToString(rankedSubsequences.get(m)));
			EMMotifFinder emmf = new EMMotifFinder(motifPositions, motifSpacing, sequences, validity);
			emmf.initializeStartingPoints(rankedSubsequences.get(m), sigmas, false, bgMean, bgVar);
			emmf.runEM(maxIterations, deltaMin);
			models.add(emmf);
			validity = eraseMotifs(validity, emmf.getStartPositions(), emmf.getMotifLength());
		}
	}
	
	private boolean checkMotifUniqueness(double[] m1, List<double[]> subsequences){

		boolean unique = true;
		if (subsequences.size() == 0) {
			return true;
		}
		
		for (double [] seq : subsequences){
			double diff = 0.0;
			for (int i=0;i<seq.length;i++){
				diff += Math.abs(seq[i] - m1[i]);
			}
			if (diff < min_motif_diff) unique = false;
			//System.out.println(getPrintableMotif(seq) + "\n" + getPrintableMotif(m1) + "\n" + diff + "\t"+ unique + "\n\n");
			
		}	
		return unique;
	}
	
	// For unit testing
	public void testUniqueCheck(){
		double [] m = {-3.318922928330818,-3.7287797352156553,-3.8228871772246875,-3.950977135300902,-4.107854052225101,-4.149026784106205};
		List<double[]> subsequences = new ArrayList<double[]>();
		
		double [] t1 = {0.0,0.0,0.0,0.0,0.0,0.0};
		double [] t2 = {-4.310223701011893,-4.701889738765324,-4.955931717267936,-5.045578197269438,-4.968974134918456,-4.707489394647542};
		double [] t3 = {0.1493956268645636,0.4145358758089242,0.29677965567829706,0.28952098298027035,0.17379503491486326,-0.17672985334843097};
		double [] t4 = {-0.9911221558409714,-0.7219020296833824,-0.7405862896140316,-0.7486392213250341,-0.8177380784618824,-0.9132493686528842};
		double [] t5 = {-4.080660337133755,-3.76324497643038,-2.9315626205232737,-1.7830852925733685,-1.3275388739723997,-1.024368751899289};
				
		subsequences.add(t1);
		subsequences.add(t2);
		subsequences.add(t3);
		subsequences.add(t4);
		subsequences.add(t5);
		
		checkMotifUniqueness(m, subsequences);

	}
	
	// For Debugging
	private String getPrintableMotif(double[] motif){
		String printable = "";
		for (double m : motif){
			printable = printable + "," + m;
		}
		
		return printable;
		
	}
	
	// make a sequence a printable string (for debugging)
	private String sequenceToString(double[] sequence){
		String printable = "";
		for (double d : sequence){
			printable = printable + d + ",";
		}
		return printable;
	}
	
	/* Get all of the subsequences of the motif length from the input sequences.
	 * Filter out ones that are too similar to eachother.
	 */
	private List<double []> getUniqueSubsequences (HashMap<String, double []> sequences, int motifPositions, int motifSpacing, int motifLength){

		System.out.println("Getting unique subsequences");
		List<double []> uniqueSubsequences = new ArrayList<double[]>();

		for (String seqName : sequences.keySet()){
			double [] seq = sequences.get(seqName);
			double [] subseq = new double [motifPositions];
			for (int i=0;i<(seq.length-motifLength);i+=(Math.ceil(motifSpacing/2.0))){
				
				// copy the candidate motif from the sequence
				for (int j=0;j<motifPositions;j++){
					subseq[j] = seq[i + j * motifSpacing];
				}
				
				//System.out.println(seqName + "\t" + getPrintableMotif(subseq)); 
				
				if (checkMotifUniqueness(subseq, uniqueSubsequences)){
					uniqueSubsequences.add(subseq);
				}
				
			}
		}
		return uniqueSubsequences;
	}
	
	/* 
	 * Run one round of EM on all subsequences and return a ranked list of starting motifs
	 */	
	private List<double []> getRankedSubsequences(	List<double []> subsequences, 
													int motifLength, int motifSpacing, double sigmaVal,
													double bgMean, double bgVar,
													HashMap<String, double []> sequences, 
													HashMap<String, boolean []> validity){
		
			System.out.println("optimizing starting motifs");
			EMMotifFinder model;
			List<double []> rankedSubsequences = new LinkedList<double []>();
			List<Double> motifProbs = new LinkedList<Double>();
			
			// initialize sigma
			double [] sigma = new double[motifLength];
			for (int s=0;s<motifLength;s++) sigma[s] = sigmaVal;
			
			// run one iteration of EM on each subsequence and calculate the likelihood
			// store them in order
			int i = 0;
			for (double[] subseq : subsequences){

				model = new EMMotifFinder(motifLength, motifSpacing, sequences, validity);
				model.initializeStartingPoints(subseq, sigma, false, bgMean, bgVar);	
				model.runEM(1, 10);
				double p = model.getLikelihoodOfMu();
				
				if (p != 0){
					int j = findInsertionIndex(motifProbs, p);
					motifProbs.add(j, p);
					rankedSubsequences.add(j, subseq);
					if (i % 100 == 0) System.out.println("tried " + i + "/" + subsequences.size());
					i++;
				}
			}

			Collections.reverse(rankedSubsequences);
			Collections.reverse(motifProbs);
			
			//for (int u=0;u<rankedSubsequences.size();u++){
			//	System.out.println(motifProbs.get(u) + "\t" + getPrintableMotif(rankedSubsequences.get(u)));
			//}
			
			return rankedSubsequences;
	}
	
	
	private int findInsertionIndex(List<Double> values, Double newValue){
		
		int max = values.size();
		int min = -1;
		int insertionIndex = min + (max - min)/2;
		
		while (max - min > 1){
			insertionIndex = min + (max - min)/2;
			if (values.get(insertionIndex) < newValue) min = insertionIndex;
			else max = insertionIndex;
		}
		
		return max;
	}
		
	/* Unimplimented */
	public void runWithSpecifiedStarts(double [] mu, int motifLength, double [] sigma, 
			HashMap<String, double []>sequences, HashMap<String, boolean []>validity){
	}
	
	public List <EMMotifFinder>  getModels(){
		return models;
	}
	
	public HashMap<String, Integer> getConsensusStartingPoints(){
		
		HashMap<String, Integer> consensusStarts = new HashMap<String, Integer>();
		
		// Count the the number of instance of each position for each iteration
		HashMap<String, HashMap<Integer, Integer>> counts = new HashMap<String, HashMap<Integer, Integer>>();
		
		for (EMMotifFinder model : models){
			Set<String> seqNames = model.getSeqNames();
			HashMap<String, Integer> startpos = model.getStartPositions();
			for (String s : seqNames){
				int sp = startpos.get(s);
				
				if (!counts.containsKey(s)){
					counts.put(s, new HashMap<Integer, Integer>());
				}
				
				HashMap<Integer, Integer> seqCounts = counts.get(s);
				
				if (!seqCounts.containsKey(sp)){
					seqCounts.put(sp, 0);
				}
				
				Integer thisCount = seqCounts.get(sp) + 1;
				seqCounts.put(sp, thisCount);
				counts.put(s, seqCounts);
			}
		}
		
		// Go through the counts and get the maxima
		for (String s : counts.keySet()){
			Integer maxCount = 0;
			Integer position = 0;
			for (Integer i : counts.get(s).keySet()){
				Integer count = counts.get(s).get(i);
				if (maxCount < count){
					position = i;
				}
			}
			consensusStarts.put(s, position);
		}
		
		return consensusStarts;
	}

	public HashMap<String, Integer> getHighestProbStartingPoints(){
		// Find the motif with the best starting probabilites.
		
		double maxP = -1;
		EMMotifFinder bestModel = models.get(0);
		
		for (EMMotifFinder model : models){
			HashMap<String, Double> spProbs = model.getStartPosProbabilities();
			Set<String> seqNames = model.getSeqNames();
			
			// multiply up all the probabilites
			double p = 1;
			for (String s : seqNames){
				p = p * (spProbs.get(s) + EPSILON);
			}
			
			if (p > maxP){
				maxP = p;
				bestModel = model;
			}
		}
		return bestModel.getStartPositions();
	}
	
	/*
	 * Use the validity track to make previously discovered motifs as invalid for the
	 * next iteration
	 */
	private HashMap<String, boolean[]> eraseMotifs(HashMap<String, boolean[]> validity,
			HashMap<String, Integer> startPositions, Integer motifLength) {
			
			HashMap<String, boolean[]> updatedValidity = new HashMap<String, boolean[]>();
		
			for (String seqName : validity.keySet()){
				
				// Copy the validity scores
				boolean[] val = validity.get(seqName);
				boolean[] updated = new boolean[val.length];
				for (int i=0;i<val.length;i++) updated[i] = val[i];
				
				// Delete the motif
				int startPos = startPositions.get(seqName);
				for (int i=0;i<motifLength-1;i++){
					updated[startPos + i] = false;
				}
				
				updatedValidity.put(seqName, updated);
			}
		
		return updatedValidity;
	}
}
