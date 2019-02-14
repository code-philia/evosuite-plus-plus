package evosuite.shell.experiment;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import evosuite.shell.FileUtils;
import evosuite.shell.excel.MergeExcels;

public class ReportTool {

	
	@Test
	public void mergeExcels() throws IOException {
		MergeExcels.excelSubfix = "_evotest_5times.xlsx";
		String reportFolder = "/Users/lylytran/Projects/Evosuite/experiments/SF100-testFilteredMethods/evoTest-reports-fbranch-14Feb";
		String outputFile = reportFolder + "/14Feb-fbranch.xlsx";
		List<String> inputFiles = FileUtils.toFilePath(MergeExcels.listExcels(reportFolder));
		MergeExcels.mergeExcel(outputFile, inputFiles, 0, false);
		
		System.out.println("Done!");
	}
	
	@Test
	public void mergeTxt() throws IOException {
		String root = "/Users/lylytran/Projects/Evosuite/modified-version/evosuite/EvosuiteTest/experiments/SF100/reports";
		List<String> inclusiveFiles = Arrays.asList(
				root + "/targetMethods-100methods.txt"
				);
		List<String> exclusivesFiles = Arrays.asList(root + "/targetMethods-invokedMethodFiltered.txt");
		
		String resultTxt = root + "/merge.txt";
		TargetMethodTool.merge(inclusiveFiles, exclusivesFiles, resultTxt);
	}
	
	@Test
	public void selectMethods() throws IOException {
		String baseDir = System.getProperty("user.dir");;
		String excelFile = baseDir + "/experiments/SF100/reports/flag-filtered-wth-GA-involved-branch.xlsx";
		String resultTxt = baseDir + "/experiments/SF100/reports/targetMethods-100methods.txt";
		String inclusiveTxt = baseDir + "/experiments/SF100/reports/targetMethods-invokedMethodFiltered.txt";
		new TargetMethodTool().selectMethods(excelFile, resultTxt, inclusiveTxt);
	}
}
