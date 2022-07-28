package org.panda.misc2.analyses;

import org.panda.resource.siteeffect.Feature;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class Zellner
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Zellner/";
	public static void main(String[] args) throws IOException
	{
		addGenesAndSites();
	}

	private static void addGenesAndSites() throws IOException
	{
		String inFile = DATA_DIR + "EC_117-M982-A01-B04-PO89-P11854-1-pSTY-PepNorm-SumInd-20210917-pep-grps-HCD-MS2-Quantitation.csv";
		String outFile = DATA_DIR + "pSTY-annotated.csv";

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "Master Protein Descriptions");
		int siteInd = ArrayUtil.indexOf(header, "Modifications in Proteins");

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("Gene\tSites\tEffect\t");
		writer.write(ArrayUtil.merge("\t", header));

		SiteEffectCollective sec = new SiteEffectCollective();

		FileUtil.linesTabbedSkip1(inFile).filter(t -> t.length > symInd).forEach(t ->
		{
			int x = t[symInd].indexOf(" GN=");
			if (x > 0)
			{
				int endInd = t[symInd].indexOf(" ", x + 4);
				String sym = endInd > 0 ? t[symInd].substring(x+4, endInd) : t[symInd].substring(x+4);

				String sites = findCanonicalSites(t[siteInd]);

				Integer effect = sec.getEffect(sym, List.of(sites.split("\\|")), Feature.PHOSPHORYLATION);

				FileUtil.lnwrite(sym + "\t" + sites + "\t" + (effect == null ? "" : effect) + "\t", writer);
			}
			else
			{
				FileUtil.lnwrite("\t\t\t", writer);
			}
			FileUtil.write(ArrayUtil.merge("\t", t), writer);
		});

		writer.close();
	}

	private static String findCanonicalSites(String positions)
	{
		if (!positions.isEmpty())
		{
			for (String s : positions.split("]; "))
			{
				String t0 = s.substring(0, s.indexOf(" "));
				if (t0.contains("Phospho")) continue;
				int endIndex = s.indexOf(" ", t0.length() + 1);
				String t1 = s.substring(t0.length()+1, endIndex);
				String t2 = s.substring(t0.length() + t1.length() + 3);

				if (!t0.contains("-") && t1.contains("Phospho"))
				{
					List<String> sites = new ArrayList<>();

					for (String tt : t2.split("; "))
					{
						if (tt.contains("(")) sites.add(tt.substring(0, tt.indexOf("(")));
					}
					return CollectionUtil.merge(sites, "|");
				}
			}
		}
		return "";
	}
}
