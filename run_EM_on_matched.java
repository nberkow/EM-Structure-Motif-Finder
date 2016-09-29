import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;


public class run_EM_on_matched {
	public static void main (String args []) throws IOException{
		
		// Read in the structure file
		StructFileReader sfr = new StructFileReader();
		
		int positions = Integer.parseInt(args[0]);
		int spacing = Integer.parseInt(args[1]);
		sfr.readFile(args[2], 39, null);
		String outPrefix = args[3];

		
		HashMap<String, double[]> scores = sfr.getStructureSequence();
		HashMap<String, boolean[]> validity = sfr.getIsValidScore();
		
		System.out.println(sfr.getMean());
		System.out.println(sfr.getVariance());
		
		EM_Manager emm = new EM_Manager();
		//emm.testUniqueCheck();
		emm.runWithOptimizedStarts(positions, spacing, 4, sfr.getVariance() * positions, 500, .00001, sfr.getVariance() * 3, sfr.getMean(), sfr.getVariance(), scores, validity);
		List<EMMotifFinder> models = emm.getModels();
		
		MotifVisualizationTool mvt = new MotifVisualizationTool();
		
		int motifNum = 1;
		for (EMMotifFinder m : models){
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(outPrefix + "_windows_motif" + motifNum + ".txt"));

			m.findMotifCenteredWinows(100);
			HashMap<String, boolean[]> valWindow = m.getValidityWindow();
			HashMap<String, double[]> scoreWindow = m.getScoreWindow();
			
			String printableWindow = mvt.getPlottableMotifWindow(scoreWindow, valWindow);
			bw.append(printableWindow);
			bw.close();
			
			bw = new BufferedWriter(new FileWriter(outPrefix + "_logo_motif" + motifNum + ".txt"));
			String printableLogo = mvt.getPlottableMotifLogo(m.mu, m.sigma);
			bw.append(printableLogo);
			bw.close();
			
			bw = new BufferedWriter(new FileWriter(outPrefix + "_startPoints_motif" + motifNum + ".txt"));
			String startPoints = m.getPrintableStartingPoints();
			bw.append(startPoints);
			bw.close();
			
			motifNum ++;
		}
	}

}
