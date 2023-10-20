package org.panda.misc2;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

public class PvsPTMPlot
{
	public static void main(String[] args)
	{
		Map<String, String> idToSym = readIDToSym();
		Set<String> ids = getIDsOfSym(idToSym, "RB1");
		ids = ids.stream().filter(id -> id.split("_").length > 2).collect(Collectors.toSet());

		Map<String, Map<String, Double>> phosMaps = readPhosphos(ids);
		phosMaps.keySet().stream().sorted().forEach(id -> System.out.println(phosMaps.get(id).size() + "\t" + id));

		Map<String, Double> pMap = readProteomic("NP_000312.2");
		Map<String, Double> ptmMap = readPhospho("NP_000312.2_S795sS807s_2_2_795_807");

		System.out.println("\nProtein\tPTM");
		Set<String> commonSamples = pMap.keySet().stream().filter(ptmMap.keySet()::contains).collect(Collectors.toSet());
		double[] x = new double[commonSamples.size()];
		double[] y = new double[commonSamples.size()];
		int i = 0;
		for (String sample : commonSamples)
		{
			x[i] = pMap.get(sample);
			y[i] = ptmMap.get(sample);
			System.out.println( x[i] + "\t" + y[i++]);
		}
		Tuple corr = Correlation.pearson(x, y);
		System.out.println("\ncorr = " + corr);
	}

	public static Map<String, String> readIDToSym()
	{
		return FileUtil.linesTabbedSkip1("/home/ozgunbabur/Data/CPTAC-PanCan/var_map_full.tsv").collect(Collectors.toMap(t -> t[0], t -> t[2]));
	}

	public static Set<String> getIDsOfSym( Map<String, String> idToSym, String sym)
	{
		return idToSym.keySet().stream().filter(id -> idToSym.get(id).equals(sym)).collect(Collectors.toSet());
	}



	public static Map<String, Double> readProteomic(String id)
	{
		return readDataFileForID("/home/ozgunbabur/Data/CPTAC-PanCan/AnalysisFiles/061721/raw/proteome_X.tsv", id);
	}

	public static Map<String, Map<String, Double>> readPhosphos(Set<String> ids)
	{
		return readDataFileForIDs("/home/ozgunbabur/Data/CPTAC-PanCan/AnalysisFiles/061721/raw/phosphoproteome_X.tsv", ids);
	}

	public static Map<String, Double> readPhospho(String id)
	{
		return readDataFileForID("/home/ozgunbabur/Data/CPTAC-PanCan/AnalysisFiles/061721/raw/phosphoproteome_X.tsv", id);
	}

	public static Map<String, Double> readDataFileForID(String file, String id)
	{
		String[] header = FileUtil.readHeader(file);
		Optional<String> first = FileUtil.lines(file).filter(l -> l.startsWith(id + "\t")).findFirst();
		if (!first.isPresent()) return null;
		String line = first.get();
		String[] t = line.split("\t");
		Map<String, Double> valMap = new HashMap<>();
		for (int i = 1; i < t.length; i++)
		{
			if (!t[i].isEmpty())
			{
				valMap.put(header[i], Double.parseDouble(t[i]));
			}
		}
		return valMap;
	}

	public static Map<String, Map<String, Double>> readDataFileForIDs(String file, Set<String> ids)
	{
		Map<String, Map<String, Double>> resultMaps = new HashMap<>();
		String[] header = FileUtil.readHeader(file);
		FileUtil.lines(file).filter(l -> ids.contains(l.substring(0, l.indexOf("\t")))).map(l -> l.split("\t")).forEach(t ->
		{
			Map<String, Double> valMap = new HashMap<>();
			for (int i = 1; i < t.length; i++)
			{
				if (!t[i].isEmpty())
				{
					valMap.put(header[i], Double.parseDouble(t[i]));
				}
			}
			resultMaps.put(t[0], valMap);
		});
		return resultMaps;
	}
}
