package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.RankedListGroupedEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class KinaseEnrichmentSignless
{
	public static void main(String[] args) throws IOException
	{
		Map<String, Map<String, Map<String, Boolean>>> rawPriors = readPriors();

//		String dataFile = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures/HRD_v2/HRDvsHRP/Full_GeneSpace/diff_expr_res/difexp/data.txt";
//		String dataFile = "/home/ozgunbabur/Analyses/Aslan-IL6/data-noadj.tsv";
//		String dataFile = "/home/ozgunbabur/Data/Aslan/gpvi_mass_spec/C2S.tsv";
		String dataFile = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/data.csv";
		List<String> rankedList = readRankedIDsFromCPFile(dataFile);
//		List<String> rankedList = readRankedIDsFromGP6Sheet(dataFile);

		// for sanity check
//		Collections.shuffle(rankedList);

		Map<String, Map<String, Set<String>>> geneToSiteToID = readDataMappingFromCPFile(dataFile);
//		Map<String, Map<String, Set<String>>> geneToSiteToID = readDataMappingFromGP6Sheet(dataFile);
		Map<String, Set<Set<String>>> priors = convertPriors(geneToSiteToID, rawPriors, rankedList);

		RankedListGroupedEnrichment.reportEnrichment(rankedList, priors, 1000000, dataFile.substring(0, dataFile.lastIndexOf(".")) + "-kinase-enrichment-v3.tsv");

		// Differential on GP6 stuff

//		String dataFile1 = "/home/ozgunbabur/Data/Aslan/gpvi_mass_spec/C1S.tsv";
//		String dataFile2 = "/home/ozgunbabur/Data/Aslan/gpvi_mass_spec/C2S.tsv";
//		List<String> rankedList1 = readRankedIDsFromGP6Sheet(dataFile1);
//		List<String> rankedList2 = readRankedIDsFromGP6Sheet(dataFile2);
//		System.out.println("rankedList1.size() = " + rankedList1.size());
//		System.out.println("rankedList2.size() = " + rankedList2.size());
//		rankedList1.retainAll(rankedList2);
//		rankedList2.retainAll(rankedList1);
//		System.out.println("After union = " + rankedList2.size());
//
//		Map<String, Map<String, Set<String>>> geneToSiteToID = readDataMappingFromGP6Sheet(dataFile1);
//		Map<String, Set<Map<String, Boolean>>> priors = convertPriors(geneToSiteToID, rawPriors, rankedList1);
//		RankedListGroupedDifferentialEnrichment.reportEnrichment(rankedList1, rankedList2, priors, 1000000,
//			dataFile1.substring(0, dataFile1.lastIndexOf(".")) + "-vs-" + dataFile2.substring(dataFile2.lastIndexOf("/") + 1, dataFile2.lastIndexOf(".")) + "-diff-kinase-enrichment.tsv");
	}

	private static Map<String, Map<String, Set<String>>> readDataMappingFromCPFile(String file)
	{
		String[] header = FileUtil.readHeader(file);
		int idInd = ArrayUtil.indexOf(header, "ID");
		int symInd = ArrayUtil.indexOf(header, "Symbols", "Genes");
		int siteInd = ArrayUtil.indexOf(header, "Sites");

		Map<String, Map<String, Set<String>>> geneToSiteToID = new HashMap<>();

		FileUtil.linesTabbedSkip1(file).filter(t -> !t[siteInd].isEmpty()).forEach(t ->
		{
			String gene = t[symInd];

			if (!geneToSiteToID.containsKey(gene)) geneToSiteToID.put(gene, new HashMap<>());

			String id = t[idInd];

			for (String site : t[siteInd].split("\\|"))
			{
				if (!geneToSiteToID.get(gene).containsKey(site)) geneToSiteToID.get(gene).put(site, new HashSet<>());
				geneToSiteToID.get(gene).get(site).add(id);
			}
		});

		return geneToSiteToID;
	}

	private static Map<String, Map<String, Set<String>>> readDataMappingFromGP6Sheet(String file)
	{
		String[] header = FileUtil.readHeader(file, 11);
		int idInd = ArrayUtil.indexOf(header, "New Sequence");
		int symInd = ArrayUtil.indexOf(header, "UniProt Gene Name");
		int siteInd = ArrayUtil.indexOf(header, "Site List");
		int uniqInd = ArrayUtil.indexOf(header, "Quan Info");

		Map<String, Map<String, Set<String>>> geneToSiteToID = new HashMap<>();

		FileUtil.linesTabbedSkip(file, 12).filter(t -> t[uniqInd].equals("Unique")).forEach(t ->
		{
			String gene = t[symInd].split(" ")[0];

			if (!geneToSiteToID.containsKey(gene)) geneToSiteToID.put(gene, new HashMap<>());

			String id = t[idInd];

			for (String site : t[siteInd].split("; "))
			{
				if (!geneToSiteToID.get(gene).containsKey(site)) geneToSiteToID.get(gene).put(site, new HashSet<>());
				geneToSiteToID.get(gene).get(site).add(id);
			}
		});

		return geneToSiteToID;
	}

//	private static String generateIDForGP6Line(String[] t, int symInd, int siteInd)
//	{
//		String gene = t[symInd].split(" ")[0];
//		String sites = t[siteInd].replaceAll("; ", "-");
//		return gene + "-" + sites;
//	}

	private static List<String> readRankedIDsFromCPFile(String file)
	{
		String[] header = FileUtil.readHeader(file);
		int idInd = ArrayUtil.indexOf(header, "ID");
		int pInd = ArrayUtil.indexOf(header, "SignedP", "HRDvsHRP_HRD", "Resting-vs-PAR4");
		int featInd = ArrayUtil.indexOf(header, "Feature", "Modification");
		Map<String, Double> idToP = FileUtil.linesTabbedSkip1(file).filter(t -> t[featInd].equals("P")).collect(Collectors.toMap(t -> t[idInd], t -> Math.abs(Double.parseDouble(t[pInd])), (d1, d2) -> d1));
		List<String> ids = new ArrayList<>(idToP.keySet());
		ids.sort(Comparator.comparing(idToP::get));
		return ids;
	}

	private static List<String> readRankedIDsFromGP6Sheet(String file)
	{
		String[] header = FileUtil.readHeader(file, 11);
		int idInd = ArrayUtil.indexOf(header, "New Sequence");
//		int symInd = ArrayUtil.indexOf(header, "UniProt Gene Name");
//		int siteInd = ArrayUtil.indexOf(header, "Site List");
//		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int uniqInd = ArrayUtil.indexOf(header, "Quan Info");

		List<String> rankedList = FileUtil.linesTabbedSkip(file, 12).filter(t -> t[uniqInd].equals("Unique")).map(t -> t[idInd]).collect(Collectors.toList());

		// remove duplicates
		Set<String> memory = new HashSet<>();
		Set<String> repeat = new HashSet<>();
		for (String s : rankedList)
		{
			if (memory.contains(s)) repeat.add(s);
			else memory.add(s);
		}
		rankedList.removeAll(repeat);

		return rankedList;
	}

	private static Map<String, Set<Set<String>>> convertPriors(Map<String, Map<String, Set<String>>> geneToSiteToID,
		Map<String, Map<String, Map<String, Boolean>>> priors, List<String> consider)
	{
		Map<String, Map<String, Set<String>>> converted = new HashMap<>();

		for (String kinase : priors.keySet())
		{
			for (String target : priors.get(kinase).keySet())
			{
				if (geneToSiteToID.containsKey(target))
				{
					for (String site : priors.get(kinase).get(target).keySet())
					{
						if (geneToSiteToID.get(target).containsKey(site))
						{
							if (!converted.containsKey(kinase)) converted.put(kinase, new HashMap<>());

							for (String id : geneToSiteToID.get(target).get(site))
							{
								if (consider.contains(id))
								{
									if (!converted.get(kinase).containsKey(target)) converted.get(kinase).put(target, new HashSet<>());
									converted.get(kinase).get(target).add(id);
								}
							}
						}
					}
				}
			}
		}

		Map<String, Set<Set<String>>> compacted = new HashMap<>();
		for (String kinase : converted.keySet())
		{
			if (converted.get(kinase).size() >= 4)
			{
				compacted.put(kinase, new HashSet<>());
				compacted.get(kinase).addAll(converted.get(kinase).values());
			}
		}

		return mergeIdenticalSets(compacted);
	}

	private static Map<String, Set<Set<String>>> mergeIdenticalSets(Map<String, Set<Set<String>>> prePriors)
	{
		Map<Set<Set<String>>, Set<String>> setsToNames = new HashMap<>();
		prePriors.forEach((name, maps) ->
		{
			if (!setsToNames.containsKey(maps)) setsToNames.put(maps, new HashSet<>());
			setsToNames.get(maps).add(name);
		});

		Map<String, Set<Set<String>>> priors = new HashMap<>();
		setsToNames.forEach((maps, names) -> priors.put(CollectionUtil.merge(names, ", "), maps));
		return priors;
	}

	private static Map<String, Map<String, Map<String, Boolean>>> readPriors()
	{
		Map<String, Map<String, Map<String, Boolean>>> priors = new HashMap<>();

		FileUtil.linesTabbed("/home/ozgunbabur/Data/causal-priors.txt")
			.filter(t -> t.length > 4 && t[1].contains("phospho") && !t[4].isEmpty()).forEach(t ->
		{
			if (!priors.containsKey(t[0])) priors.put(t[0], new HashMap<>());
			Map<String, Map<String, Boolean>> groups = priors.get(t[0]);
			if (!groups.containsKey(t[2])) groups.put(t[2], new HashMap<>());
			Map<String, Boolean> group = groups.get(t[2]);

			boolean direction = t[1].startsWith("p");

			for (String site : t[4].split(";"))
			{
				group.put(site, direction);
			}
		});
		return priors;
	}
}
