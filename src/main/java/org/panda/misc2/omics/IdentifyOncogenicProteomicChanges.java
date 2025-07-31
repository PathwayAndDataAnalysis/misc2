package org.panda.misc2.omics;

import org.panda.misc2.causalpath.CausalPathSubnetwork;
import org.panda.resource.GeminiLungCancerGenes;
import org.panda.resource.OncoKB;
import org.panda.resource.PCPathway;
import org.panda.resource.PCPathwayHGNC;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.siteeffect.Feature;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.TTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * In a comparison-based setting, using t-test
 */
public class IdentifyOncogenicProteomicChanges
{
	private static SiteEffectCollective sec = new SiteEffectCollective();

	public static void main(String[] args) throws IOException
	{
		exploreCPTACLSCCClustersVsNormals();
	}

	public static void test() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/";
		Map<String[], Tuple> changedRows = getChangedRows(
			dir + "data.txt", dir + "tumor-vs-normal/parameters.txt", 0.1);
		Map<String[], Tuple> oncogenicRows = selectOncogenicRows(changedRows);
		System.out.println("oncogenicRows.size() = " + oncogenicRows.size());
		Set<String> genes = oncogenicRows.keySet().stream().map(r -> r[1].split(" ")).flatMap(Arrays::stream).collect(Collectors.toSet());
		System.out.println("genes.size = " + genes.size());
		oncogenicRows.keySet().forEach(r -> System.out.println(Arrays.toString(r)));
	}

	private static void exploreCPTACLSCCClustersVsNormals() throws IOException
	{
		String base = "/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/";
		String dir = base + "clusters/against-normals/";

		List<String> names = new ArrayList<>();
		List<Set<String>> idSets = new ArrayList<>();
		Set<String> allGenes = new HashSet<>();
		Set<String> focusGenes = new HashSet<>();

		// Remember each oncogenic change for each subtype
		Map<String, Map<String[], Tuple>> tupleMap = new HashMap<>();

		for (File d : new File(dir).listFiles())
		{
			if (d.toString().endsWith(".DS_Store")) continue;

			String name = d.toString().substring(d.toString().lastIndexOf("/") + 1);
			Map<String[], Tuple> changedRows = getChangedRows(
				base + "data.txt", d.toString() + "/parameters.txt", 0.1);
			Map<String[], Tuple> oncogenicRows = selectOncogenicRows(changedRows);
			tupleMap.put(d.toString().substring(d.toString().lastIndexOf("/") + 1), oncogenicRows);
			oncogenicRows.keySet().stream().map(s -> s[1]).forEach(g -> focusGenes.addAll(Arrays.asList(g.split(" "))));

			Set<String> ids = oncogenicRows.keySet().stream().map(k -> k[0]).collect(Collectors.toSet());

			names.add(name);
			idSets.add(ids);

			Map<String[], Tuple> allRows = readRows(base + "data.txt", d.toString() + "/parameters.txt");
			allRows.keySet().stream().map(t -> t[1]).forEach(g -> allGenes.addAll(Arrays.asList(g.split(" "))));
		}

		CollectionUtil.printNameMapping(names.toArray(new String[0]));
		CollectionUtil.printVennSets(idSets.toArray(new Collection[0]));

		Set<String> all = CollectionUtil.getUnion(idSets.toArray(new Collection[0]));
		System.out.println("all.size() = " + all.size() + "\n");

		printAsList(tupleMap);
		printEnrichment(focusGenes, allGenes);
	}

	private static Map<String[], Tuple> selectOncogenicRows(Map<String[], Tuple> rows)
	{
		return rows.keySet().stream().filter(r -> isChangeOncogenic(r, rows.get(r))).collect(Collectors.toMap(r -> r, rows::get));
	}

	private static Map<String[], Tuple> getChangedRows(String matrixFile, String parametersFile, double fdrThr) throws IOException
	{
		Map<String[], Tuple> ttMap = readRows(matrixFile, parametersFile);

//		Map<String[], Double> pvals = ttMap.keySet().stream().collect(Collectors.toMap(k -> k, k -> ttMap.get(k).p));
		Map<String[], Double> pvalsP = ttMap.keySet().stream().filter(k -> !k[2].isEmpty()).collect(Collectors.toMap(k -> k, k -> ttMap.get(k).p));
		Map<String[], Double> pvalsT = ttMap.keySet().stream().filter(k -> k[2].isEmpty()).collect(Collectors.toMap(k -> k, k -> ttMap.get(k).p));

//		List<String[]> select = FDR.select(pvals, null, fdrThr);
		List<String[]> selectP = FDR.select(pvalsP, null, fdrThr);
		List<String[]> selectT = FDR.select(pvalsT, null, fdrThr);
		List<String[]> select = new ArrayList<>(selectP);
		select.addAll(selectT);

		return select.stream().collect(Collectors.toMap(i -> i, ttMap::get));
	}

	private static Map<String[], Tuple> readRows(String matrixFile, String parametersFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(matrixFile)).findFirst().get().split("\t");
		List<String> testSamples = readTestSamples(parametersFile);
		List<String> ctrlSamples = readControlSamples(parametersFile);
		int[] tInd = getIndexes(header, testSamples);
		int[] cInd = getIndexes(header, ctrlSamples);

		Map<String[], Tuple> ttMap = new HashMap<>();

		Files.lines(Paths.get(matrixFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String[] id = new String[4];
			System.arraycopy(t, 0, id, 0, id.length);
			double[] test = extractVals(t, tInd);
			double[] ctrl = extractVals(t, cInd);

			Tuple tuple = TTest.testPaired(ctrl, test);
			if (!tuple.isNaN()) ttMap.put(id, tuple);
		});
		return ttMap;
	}


	private static List<String> readTestSamples(String parametersFile) throws IOException
	{
		return Files.lines(Paths.get(parametersFile)).filter(l -> l.startsWith("test-value-column = "))
			.map(l -> l.substring(l.indexOf("=") + 1).trim()).collect(Collectors.toList());
	}

	private static List<String> readControlSamples(String parametersFile) throws IOException
	{
		return Files.lines(Paths.get(parametersFile)).filter(l -> l.startsWith("control-value-column = "))
			.map(l -> l.substring(l.indexOf("=") + 1).trim()).collect(Collectors.toList());
	}

	private static int[] getIndexes(String[] header, List<String> cols)
	{
		int[] ind = new int[cols.size()];
		for (int i = 0; i < ind.length; i++)
		{
			ind[i] = ArrayUtil.indexOf(header, cols.get(i));
		}
		return ind;
	}

	private static double[] extractVals(String[] row, int[] inds)
	{
		double[] vals = new double[inds.length];
		for (int i = 0; i < vals.length; i++)
		{
			vals[i] = Double.valueOf(row[inds[i]]);
		}
		return vals;
	}

	private static Integer getOncogeneOrTumorSupp(String[] row)
	{
		String[] genes = row[1].split(" ");

		for (int i = 0; i < genes.length; i++)
		{
			if (OncoKB.get().isOncogeneOnly(genes[i])) return 1;
			else if (OncoKB.get().isTumorSuppressorOnly(genes[i])) return -1;
		}
		return null;
	}

	private static Integer getActivatingOrInhibiting(String[] row)
	{
		String[] genes = row[1].split(" ");
		for (int i = 0; i < genes.length; i++)
		{
			if (!row[2].isEmpty())
			{
				for (String site : row[2].split(" ")[i].split("\\|"))
				{
					Integer effect = sec.getEffect(genes[i], site, Feature.PHOSPHORYLATION);
					if (effect != null && effect != 0)
					{
						return effect;
					}
				}
			}
		}
		return null;
	}

	private static boolean isChangeOncogenic(String[] row, Tuple tup)
	{
		String[] genes = row[1].split(" ");
		boolean up = tup.v > 0;
		for (int i = 0; i < genes.length; i++)
		{
			if (row[2].isEmpty())
			{
				if (isChangeOncogenic(genes[i], null, null, up)) return true;
			}
			else
			{
				for (String site : row[2].split(" ")[i].split("\\|"))
				{
					if (isChangeOncogenic(genes[i], site, Feature.getFeat(row[3]), up)) return true;
				}
			}
		}
		return false;
	}

	public static boolean isChangeOncogenic(String gene, String site, Feature mod, boolean up)
	{
		if (OncoKB.get().isCancerGene(gene))
		{
			if (site == null)
			{
				if (up) return OncoKB.get().isOncogeneOnly(gene);
				else return OncoKB.get().isTumorSuppressorOnly(gene);
			}
			else
			{
				Integer effect = sec.getEffect(gene, site, mod);
				if (effect != null && effect != 0)
				{
					if (effect == 1)
					{
						if (up) return OncoKB.get().isOncogeneOnly(gene);
						else return OncoKB.get().isTumorSuppressorOnly(gene);
					}
					else if (effect == -1)
					{
						if (up) return OncoKB.get().isTumorSuppressorOnly(gene);
						else return OncoKB.get().isOncogeneOnly(gene);
					}
				}
			}
		}

		return false;
	}

	public static boolean isChangeOncogenicForLung(String gene, String site, Feature mod, boolean up)
	{
		if (GeminiLungCancerGenes.get().isCancerGene(gene))
		{
			if (site == null)
			{
				if (up) return GeminiLungCancerGenes.get().isOncogeneOnly(gene);
				else return GeminiLungCancerGenes.get().isTumorSuppressorOnly(gene);
			}
			else
			{
				Integer effect = sec.getEffect(gene, site, mod);
				if (effect != null && effect != 0)
				{
					if (effect == 1)
					{
						if (up) return GeminiLungCancerGenes.get().isOncogeneOnly(gene);
						else return GeminiLungCancerGenes.get().isTumorSuppressorOnly(gene);
					}
					else if (effect == -1)
					{
						if (up) return GeminiLungCancerGenes.get().isTumorSuppressorOnly(gene);
						else return GeminiLungCancerGenes.get().isOncogeneOnly(gene);
					}
				}
			}
		}

		return false;
	}

	// select oncogenic AND tarrgetable, otherwise use the above one.
//	public static boolean isChangeOncogenic(String gene, String site, boolean up)
//	{
//		if (OncoKB.get().isCancerGene(gene))
//		{
//			if (site == null)
//			{
//				if (up) return OncoKB.get().isOncogeneOnly(gene);
//			}
//			else
//			{
//				Integer effect = sec.getEffect(gene, site);
//				if (effect != null && effect != 0)
//				{
//					if (effect == 1)
//					{
//						if (up) return OncoKB.get().isOncogeneOnly(gene);
//					}
//					else if (effect == -1)
//					{
//						if (!up) return OncoKB.get().isOncogeneOnly(gene);
//					}
//				}
//			}
//		}
//
//		return false;
//	}

	public static Integer getSiteEffect(String gene, String site)
	{
		return sec.getEffect(gene, site, Feature.PHOSPHORYLATION);
	}

	private static void printAsList(Map<String, Map<String[], Tuple>> tupleMap)
	{
		Map<String, Set<String>> idToSubtype = new HashMap<>();
		Map<String, String[]> idToRow = new HashMap<>();
		Map<String, Tuple> idToTup = new HashMap<>();

		List<String> subtypes = new ArrayList<>();

		for (String subtype : tupleMap.keySet())
		{
			Map<String[], Tuple> tups = tupleMap.get(subtype);
			for (String[] quad : tups.keySet())
			{
				String id = quad[0];
				if (!idToSubtype.containsKey(id)) idToSubtype.put(id, new HashSet<>());
				idToSubtype.get(id).add(subtype);
				if (!subtypes.contains(subtype)) subtypes.add(subtype);

				Tuple tup = tups.get(quad);

				// remember the most significant p
				if (!idToTup.containsKey(id) || idToTup.get(id).p > tup.p) idToTup.put(id, tup);

				idToRow.put(id, quad);
			}
		}

		System.out.print("ID\tGene Type\tMeasurement type\tChange");
		subtypes.sort(String::compareTo);
		subtypes.forEach(st -> System.out.print("\t" + st));
		System.out.print("\tMost sig. P val.");

		idToRow.keySet().stream().sorted((k1, k2) ->
		{
//			int c = Integer.compare(idToSubtype.get(k2).size(), idToSubtype.get(k1).size());
			int c = Double.compare(idToTup.get(k1).p, idToTup.get(k2).p);
			if (c == 0)
			{
				c = k1.compareTo(k2);
			}
			return c;
		}).forEach(id ->
		{
			System.out.print("\n" + id);
			Tuple tup = idToTup.get(id);
			String[] row = idToRow.get(id);
			Integer onco = getOncogeneOrTumorSupp(row);
			Integer pType = getActivatingOrInhibiting(row);
			System.out.print("\t" + (onco==null ? "" : onco==1 ? "Oncogene" : onco==-1 ? "Tumor Suppressor" : ""));
			System.out.print("\t" + (pType==null ? "Global protein" : pType==1 ? "Activating phosphosite" : pType==-1 ? "Inhibiting phosphosite" : "Error"));
			System.out.print("\t" + (tup.v > 0 ? "Upregulated" : "Downregulated"));

			subtypes.forEach(st -> System.out.print("\t" + (idToSubtype.get(id).contains(st) ? "X" : "")));
			System.out.print("\t" + idToTup.get(id).p);
		});
	}

	private static void printEnrichment(Set<String> focusGenes, Set<String> allGenes) throws IOException
	{
		PCPathway pc = new PCPathway();
		pc.writeEnrichmentResults(focusGenes, allGenes, 5, 200, "/Users/ozgun/Documents/Temp/enrich.txt");

	}
}
