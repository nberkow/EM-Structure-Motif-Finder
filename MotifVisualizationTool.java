import java.util.HashMap;


public class MotifVisualizationTool {
	
	public String getPlottableMotifWindow(
				HashMap<String, double[]> motifCenteredScoreWindow, 
				HashMap<String, boolean[]> motifCenteredValidityWindow){
		
		String plottableData = "";
		
		for (String s : motifCenteredScoreWindow.keySet()){
			String reportLine = s;
			double [] scores = motifCenteredScoreWindow.get(s);
			boolean [] validity = motifCenteredValidityWindow.get(s);
			for (int i=0;i<scores.length;i++){
				if (validity[i]) {
					reportLine = reportLine + "," + scores[i];
				}
				else {
					reportLine = reportLine + ",NA";
				}
			}
			plottableData = plottableData + reportLine + "\n";
		}
		
		return plottableData;
	}
	
	public String getPlottableMotifLogo(double[] mu, double [] sigma){
		String logo = "";
		for (double m : mu){
			logo = logo + "," + m;
		}
		logo = logo + "\n";
		for (double s : sigma){
			logo = logo + "," + s;
		}
		return logo;
	}
}
