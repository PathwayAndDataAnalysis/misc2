package org.panda.misc2.analyses;

import org.panda.utility.FileUtil;

import java.util.HashMap;
import java.util.Map;

public class Ghatge
{
	public static final String BASE = "/home/ozgunbabur/Analyses/Platelet-Blood-paper/ghatge-metabolome/";
	public static void main(String[] args)
	{
		readChEBIMap();
	}

	public static Map<String, Double> readChEBIMap()
	{
		Map<String, String> nameToChEBI = new HashMap<>();
		FileUtil.lines(BASE + "query-with-name.csv").skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0].replaceAll("\"", "");
			String cheBI = t[4].replaceAll("\"", "");
			if (!cheBI.equals("NA"))
			{
				nameToChEBI.put(name, cheBI);
			}
		});
		Map<String, String> pubchemToChEBI = new HashMap<>();
		FileUtil.lines(BASE + "query-with-pubchem.csv").skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String pubchem = t[0].replaceAll("\"", "");
			String cheBI = t[4].replaceAll("\"", "");
			if (!cheBI.equals("NA"))
			{
				pubchemToChEBI.put(pubchem, cheBI);
			}
		});

		Map<String, Double> chebiToFC = new HashMap<>();
		FileUtil.linesTabbedSkip1(BASE + "Supp-Table1.csv").filter(t -> t.length > 5).forEach(t ->
		{
			String name = t[0];
			String pubchem = t[3];

			String chebi = null;

			if (!pubchem.isEmpty()) chebi = pubchemToChEBI.get(pubchem);
			if (chebi == null) chebi = nameToChEBI.get(name);

			if (chebi != null && !t[5].isEmpty())
			{
				double fc = Double.parseDouble(t[4]);
				double q = Double.parseDouble(t[5]);

				if (q > 0.1) fc = 0;
				else
				{
					if (fc >= 1) fc = fc - 1;
					else
					{
						fc = -((1 / fc) - 1);
					}
				}

				chebiToFC.put("CHEBI:" + chebi, fc);
			}
		});

		chebiToFC.forEach((id, fc) -> System.out.println(id + "\t" + id + "\t\t\tC\t" + fc));

		return chebiToFC;
	}
}
