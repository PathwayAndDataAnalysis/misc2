package org.panda.misc2.causalpath;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This class is for CausalPath-formatted data that has multiple SignedP columns. It assumes difference in significance
 * is an indicator for difference in amounts. This may be the case when the tests are of similar quality, with similar
 * amounts of noise, and similar statistical strength.
 */
public class CausalPathDifferenceNetwork
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/";
		List<String> newNames = new ArrayList<>(Arrays.asList("PAR4-vs-PAR1", "Thrombin-vs-PAR4", "Thrombin-vs-PAR1"));
		List<String[]> comparisons = new ArrayList<>(Arrays.asList(
			new String[]{"Resting-vs-PAR4", "Resting-vs-PAR1"},
			new String[]{"Resting-vs-Thrombin", "Resting-vs-PAR4"},
			new String[]{"Resting-vs-Thrombin", "Resting-vs-PAR1"}));
		generateDiffData(dir + "data.csv", dir + "data-diff.csv", newNames, comparisons, 5);
	}

	public static void generateDiffData(String inFile, String outFile, List<String> newNames, List<String[]> comparisons,
		int numberOfPreservedColumns) throws IOException
	{
		String[] header = FileUtil.readHeader(inFile);
		List<int[]> comparisonInds = new ArrayList<>();
		for (int i = 0; i < comparisons.size(); i++)
		{
			int[] inds = new int[2];
			inds[0] = ArrayUtil.indexOf(header, comparisons.get(i)[0]);
			inds[1] = ArrayUtil.indexOf(header, comparisons.get(i)[1]);
			comparisonInds.add(inds);
		}
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		for (int i = 0; i < numberOfPreservedColumns; i++)
		{
			if (i > 0) FileUtil.write("\t", writer);
			writer.write(header[i]);
		}
		for (String newName : newNames)
		{
			writer.write("\t" + newName);
		}

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			FileUtil.write("\n", writer);

			for (int i = 0; i < numberOfPreservedColumns; i++)
			{
				if (i > 0) FileUtil.write("\t", writer);
				FileUtil.write(t[i], writer);
			}
			for (int[] ind : comparisonInds)
			{
				int sign0 = t[ind[0]].startsWith("-") ? -1 : 1;
				int sign1 = t[ind[1]].startsWith("-") ? -1 : 1;
				double v0 = Math.abs(Double.parseDouble(t[ind[0]]));
				double v1 = Math.abs(Double.parseDouble(t[ind[1]]));
				v0 = -Math.log10(v0) * sign0;
				v1 = -Math.log10(v1) * sign1;
				double diff = v0 - v1;
				FileUtil.tab_write(diff, writer);
			}
		});

		writer.close();
	}
}
