package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

public class EnrichmentResultDependencyGrapher
{
	public static void main(String[] args) throws IOException
	{
//		String[] cases = new String[]{"male-vs-female", "old-female-vs-young-female", "old-male-vs-old-female", "old-male-vs-young-male", "old-vs-young", "young-male-vs-young-female"};
//
//		for (String aCase : cases)
//		{
//			String dir = "/home/ozgunbabur/Analyses/Aslan-YoungAndOld/WithKinaseLib/" + aCase + "/";
//			generateJaccardDistGraph(dir + "kinase-enrichment-specific.tsv", "Name", "Direction",
//				"P-value", "FDR", "Members in top 25%", 0.1, 0.6,
//				dir + "dependency-graph-specific");
//			generateJaccardDistGraph(dir + "kinase-enrichment.tsv", "Name", "Direction",
//				"P-value", "FDR", "Members in top 25%", 0.1, 0.6,
//				dir + "dependency-graph");
//		}

		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/WithKinaseLibrary/Resting-vs-Thrombin/";
		generateJaccardDistGraph(dir + "kinase-enrichment-specific.tsv", "Name", "Direction",
			"P-value", "FDR", "Members in top 25%", 0.1, 0.6,
			dir + "dependency-graph-specific");


	}

	public static void generateJaccardDistGraph(String inputFile, String nameCol, String signCol, String pCol,
												String fdrCol, String membersCol, double fdrThr, double jaccardThr,
												String outSIFWOExt) throws IOException
	{
		System.out.println("inputFile = " + inputFile);
		String[] header = FileUtil.readHeader(inputFile);
		int nameInd = ArrayUtil.indexOf(header, nameCol);
		int signInd = ArrayUtil.indexOf(header, signCol);
		int pInd = ArrayUtil.indexOf(header, pCol);
		int fdrInd = ArrayUtil.indexOf(header, fdrCol);
		int membersInd = ArrayUtil.indexOf(header, membersCol);

		Map<String, Integer> signMap = new HashMap<>();
		Map<String, Double> pMap = new HashMap<>();
		Map<String, Double> fdrMap = new HashMap<>();
		Map<String, Set<String>> membersMap = new HashMap<>();

		FileUtil.linesTabbedSkip1(inputFile).forEach(t ->
		{
			String name = t[nameInd];
			signMap.put(name, t[signInd].equals("+") ? 1 : -1);
			pMap.put(name, Double.parseDouble(t[pInd]));
			fdrMap.put(name, Double.parseDouble(t[fdrInd]));
			membersMap.put(name, t.length > membersInd ? new HashSet<>(Arrays.asList(t[membersInd].split(" "))) : Collections.emptySet());
		});

		double pThr = 0;
		for (String name : fdrMap.keySet())
		{
			double fdr = fdrMap.get(name);
			double p = pMap.get(name);
			if (fdr <= fdrThr & p > pThr)
			{
				pThr = p;
			}
		}

		BufferedWriter sifWriter = FileUtil.newBufferedWriter(outSIFWOExt + ".sif");
		BufferedWriter fmtWriter = FileUtil.newBufferedWriter(outSIFWOExt + ".format");
		fmtWriter.write("node\tall-nodes\tbordercolor\t100 100 100\n" +
			"edge\tall-edges\twidth\t2");

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{new Color(100, 100, 250), Color.WHITE, new Color(250, 100, 100)});

		for (String name : pMap.keySet())
		{
			double pval = pMap.get(name);
			if (pval <= pThr)
			{
				double val = -Math.log(pval) * signMap.get(name);
				FileUtil.lnwrite("node\t" + name + "\tcolor\t" + vtc.getColorInString(val), fmtWriter);
			}
		}

		vtc = new ValToColor(new double[]{jaccardThr - 0.1, 1}, new Color[]{Color.WHITE, Color.BLACK});
		Set<String> connected = new HashSet<>();

		for (String n1 : pMap.keySet())
		{
			if (pMap.get(n1) <= pThr)
			{
				for (String n2 : pMap.keySet())
				{
					if (n1.compareTo(n2) < 0 && pMap.get(n2) < pThr)
					{
						Set<String> m1 = membersMap.get(n1);
						Set<String> m2 = membersMap.get(n2);
						double jac = CollectionUtil.getJaccardSimilarity(m1, m2);

						if (jac >= jaccardThr)
						{
							FileUtil.writeln(n1 + "\tinteracts-with\t" + n2, sifWriter);
							connected.add(n1);
							connected.add(n2);
							FileUtil.lnwrite("edge\t" + n1 + " interacts-with " + n2 + "\tcolor\t" + vtc.getColorInString(jac), fmtWriter);
						}
					}
				}
			}
		}

		for (String name : pMap.keySet())
		{
			if (!connected.contains(name) && pMap.get(name) <= pThr) FileUtil.writeln(name, sifWriter);
		}

		sifWriter.close();
		fmtWriter.close();
	}
}
