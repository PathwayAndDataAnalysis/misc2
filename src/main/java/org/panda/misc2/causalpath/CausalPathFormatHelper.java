package org.panda.misc2.causalpath;

import org.panda.resource.OncoKB;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.siteeffect.Feature;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.util.*;
import java.util.List;

public class CausalPathFormatHelper
{
	public static void main(String[] args)
	{
//		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-PanCan/";
//
//		List<String> lines = getFormatLinesAdjSignedP(new HashSet<>(Arrays.asList("CDK1", "CDK2", "RB1")), dir + "prot-data-diffexp.txt",
//			dir + "rna-data-diffexp.txt", "ID", "Symbols", "Sites", "Modification",
//			"Effect", "24", 0.1);
//
//		lines.forEach(System.out::println);

		String dir = "/Users/ozgun/Documents/Analyses/WilliamGrey/";
		extractSymbols(dir + "causative.format", dir + "symbols.txt");
	}

	public static List<String> getFormatLinesAdjSignedP(Set<String> genes, String proteinData, String rnaData,
		String idCol, String genesCol, String sitesCol, String modCol, String effectCol, String pCol, double sigThr)
	{
		String[] pHeader = FileUtil.readHeader(proteinData);
		String[] rHeader = FileUtil.readHeader(rnaData);

		return getFormatLines(genes, proteinData, rnaData,
			row -> row[ArrayUtil.indexOf(pHeader, idCol)],
			row -> row[ArrayUtil.indexOf(pHeader, genesCol)],
			row -> row[ArrayUtil.indexOf(pHeader, sitesCol)],
			row -> row[ArrayUtil.indexOf(pHeader, modCol)],
			row -> row[ArrayUtil.indexOf(pHeader, effectCol)],
			new SignedPValSelect(ArrayUtil.indexOf(pHeader, pCol), sigThr),
			new SignedPValSelect(ArrayUtil.indexOf(rHeader, pCol), sigThr),
			new SignedPValChange(ArrayUtil.indexOf(pHeader, pCol)),
			new SignedPValChange(ArrayUtil.indexOf(rHeader, pCol)));
	}

	public static List<String> getFormatLines(Set<String> genes, String proteinData, String rnaData,
		CellReader idReader, CellReader genesReader, CellReader sitesReader, CellReader modReader, CellReader effectReader,
		Selector pSelector, Selector rSelector, Change pChange, Change rChange)
	{
		List<String> lines = new ArrayList<>();

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10},
			new Color[]{new Color(40, 80, 255), Color.WHITE, new Color(255, 80, 40)});

		SiteEffectCollective sec = new SiteEffectCollective();

		FileUtil.lines(proteinData).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String gene = genesReader.read(t);

			if (genes.contains(gene) && pSelector.select(t))
			{
				if (modReader.read(t).isEmpty())
				{
					lines.add("node\t" + gene + "\tcolor\t" + vtc.getColorInString(pChange.getChangeValue(t)));
				}
				else
				{
					String lett = modReader.read(t).toLowerCase();

					String e = effectReader.read(t);
					Integer effect = e.isEmpty() || e.equals("c") ? null : e.equals("a") ? 1 : -1;

					if (effect == null)
					{
						effect = sec.getEffect(gene, Arrays.asList(sitesReader.read(t).split(";")),
							Feature.getFeat(modReader.read(t)));
					}

					String bordColor = effect == null || effect == 0 ? "0 0 0" : effect == 1 ? "0 180 20" : "180 0 20";

					lines.add("node\t" + gene + "\trppasite\t" + idReader.read(t) + "|" + lett + "|" +
						vtc.getColorInString(pChange.getChangeValue(t)) + "|" + bordColor + "|" +
						pChange.getChangeValue(t));
				}
			}
		});

//		FileUtil.lines(rnaData).skip(1).map(l -> l.split("\t")).forEach(t ->
//		{
//			if (genes.contains(t[0]) && rSelector.select(t))
//			{
//				lines.add("node\t" + t[0] + "\trppasite\t" + t[0] + "-rna|r|" +
//					vtc.getColorInString(rChange.getChangeValue(t)) + "|0 0 0|" +
//					rChange.getChangeValue(t));
//			}
//		});

		return lines;
	}

	public static Set<String> getFormatLinesForOncogenicEffect(Set<String> genes)
	{
		Set<String> lines = new HashSet<>();

		for (String gene : genes)
		{
			if (OncoKB.get().isCancerGene(gene))
			{
				if (OncoKB.get().isOncogeneOnly(gene))
				{
					lines.add("node\t" + gene + "\tbordercolor\t255 50 50");
				}
				else if (OncoKB.get().isTumorSuppressorOnly(gene))
				{
					lines.add("node\t" + gene + "\tbordercolor\t50 50 255");
				}
				else
				{
					lines.add("node\t" + gene + "\tbordercolor\t180 180 50");
				}
			}
		}

		return lines;
	}

	public interface CellReader
	{
		String read(String[] row);
	}

	public interface Selector
	{
		boolean select(String[] row);
	}

	public interface Change
	{
		double getChangeValue(String[] row);
	}

	/**
	 * A change detector for adjusted signed p-values
	 */
	static class SignedPValChange implements Change
	{
		int pIndex;

		public SignedPValChange(int pIndex)
		{
			this.pIndex = pIndex;
		}

		@Override
		public double getChangeValue(String[] row)
		{
			double v = row.length <= pIndex || row[pIndex].isEmpty() || row[pIndex].equals("NaN") || row[pIndex].equals("NA") ? Double.NaN : Double.valueOf(row[pIndex]);
			int s = (int) Math.signum(v);
			v = Math.abs(v);
			v = -Math.log(v);
			v *= s;
			return v;
		}
	}
	/**
	 * A selector for adjusted signed p-values
	 */
	static class SignedPValSelect implements Selector
	{
		int pIndex;
		double thr;

		public SignedPValSelect(int pIndex, double thr)
		{
			if (pIndex < 0)
			{
				System.out.println();
			}
			this.pIndex = pIndex;
			this.thr = thr;
		}

		@Override
		public boolean select(String[] row)
		{
			double v = row.length <= pIndex || row[pIndex].isEmpty() || row[pIndex].equals("NaN") || row[pIndex].equals("NA") ? Double.NaN : Double.valueOf(row[pIndex]);
			v = Math.abs(v);
			return v <= thr;
		}
	}

	static void extractSymbols(String formatFile, String outSymbolFile)
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(outSymbolFile);

		FileUtil.linesTabbed(formatFile).filter(t -> t[0].equals("node") && !t[1].equals("all-nodes")).map(t -> t[1])
			.distinct().forEach(sym -> FileUtil.writeln(sym, writer));

		FileUtil.closeWriter(writer);
	}
}
