package org.panda.misc2.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PNNLTCell
{
	public static final String DIR = "/home/ozgunbabur/Analyses/PNNLTCell/";

	public static void main(String[] args)
	{
		mapToHuman(DIR + "mouseTcell20min.csv", DIR + "20-min/humanized-20min.tsv");
	}

	public static void mapToHuman(String mFile, String hFile)
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(hFile);
		FileUtil.write("ID\tSymbols\tSites\tEffect\tFeature\tSignedP", writer);

		String[] header = FileUtil.readHeader(mFile);
		int upIndex = ArrayUtil.indexOf(header, "T: ProteinID");
		int siteIndex = ArrayUtil.indexOf(header, "STY sites");
		int pIndex = ArrayUtil.indexOf(header, "N: -Log ANOVA p value");
		int dirIndex = ArrayUtil.indexOf(header, "Direction");

		FileUtil.linesTabbedSkip1(mFile).forEach(t ->
		{
			String mUP = t[upIndex];
			String[] mSites = parseMSites(t[siteIndex]);

			Map<String, List<String>> mappedSites = SiteMappingMouseToHuman.get().mapToHumanSite(mUP, mSites);

			if (!mappedSites.isEmpty())
			{
				List<String> syms = new ArrayList<>();
				List<String> sites = new ArrayList<>();
				for (String hUP : mappedSites.keySet())
				{
					String sym = HGNC.get().getSymbol(hUP);
					if (sym != null)
					{
						syms.add(sym);
						sites.add(CollectionUtil.merge(mappedSites.get(hUP), "|"));
					}
				}

				String symbolsStr = CollectionUtil.merge(syms, " ");
				String sitesStr = CollectionUtil.merge(sites, " ");

				double p = Math.pow(10, -Double.parseDouble(t[pIndex]));
				if (t[dirIndex].startsWith("D")) p = -p;

				String id = generateID(symbolsStr, sitesStr);
				FileUtil.lnwrite(id + "\t" + symbolsStr + "\t" + sitesStr + "\t\tP\t" + p, writer);
			}

		});

		FileUtil.closeWriter(writer);
	}

	static String generateID(String symbols, String sites)
	{
		String[] sy = symbols.split(" ");
		String[] si = sites.split(" ");
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < sy.length; i++)
		{
			sb.append(sy[i]).append("-").append(si[i].replaceAll("\\|", "-")).append("-");
		}
		String result = sb.toString();
		result = result.substring(0, result.length() - 1);
		return result;
	}
	public static String[] parseMSites(String siteStr)
	{
		List<Integer> letInds = new ArrayList<>();
		for (int i = 0; i < siteStr.length(); i++)
		{
			if(!Character.isDigit(siteStr.charAt(i))) letInds.add(i);
		}
		String[] sites = new String[letInds.size()];

		for (int i = 0; i < letInds.size(); i++)
		{
			String site = i == letInds.size() - 1 ?
				siteStr.substring(letInds.get(i)) :
				siteStr.substring(letInds.get(i), letInds.get(i + 1));

			sites[i] = site;
		}

		return sites;
	}
}
