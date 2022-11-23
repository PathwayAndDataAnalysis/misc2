package org.panda.misc2;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class TestDataSingleCell
{
	static final double LOW = 2;
	static final double HIGH = 10;
	static final double SD = 2;
	static final double SD_BRIDGE = 2;
	static final int C1_SIZE = 100;
	static final int C2_SIZE = 100;
	static final int C3_SIZE = 100;
	static final int C1_C2_BRIDGE_SIZE = 50;

	static final int P1_GENE_SIZE = 20;
	static final int P2_GENE_SIZE = 20;
	static final int P3_GENE_SIZE = 20;

	static Random rand = new Random();

	public static void main(String[] args) throws IOException
	{
		generate("/home/ozgunbabur/Downloads/temp/test-data.tsv");
	}

	public static void generate(String outFile) throws IOException
	{
		List<String> cells = getCellNames();
		List<String> genes = getGeneNames();

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);

		writer.write("Gene");
		cells.forEach(c -> FileUtil.tab_write(c, writer));

		for (String gene : genes)
		{
			FileUtil.lnwrite(gene, writer);

			for (String cell : cells)
			{
				FileUtil.tab_write(generateValueFor(gene, cell), writer);
			}
		}

		writer.close();
	}

	static List<String> getGeneNames()
	{
		List<String> genes = new ArrayList<>();
		for (int i = 0; i < P1_GENE_SIZE; i++)
		{
			genes.add("G1_" + i);
		}
		for (int i = 0; i < P2_GENE_SIZE; i++)
		{
			genes.add("G2_" + i);
		}
		for (int i = 0; i < P3_GENE_SIZE; i++)
		{
			genes.add("G3_" + i);
		}
		return genes;
	}

	static List<String> getCellNames()
	{
		List<String> cells = new ArrayList<>();
		for (int i = 0; i < C1_SIZE; i++)
		{
			cells.add("C1_" + i);
		}
		for (int i = 0; i < C2_SIZE; i++)
		{
			cells.add("C2_" + i);
		}
		for (int i = 0; i < C3_SIZE; i++)
		{
			cells.add("C3_" + i);
		}
		for (int i = 0; i < C1_C2_BRIDGE_SIZE; i++)
		{
			cells.add("B_" + i);
		}
		return cells;
	}

	static double generateValueFor(String gene, String cell)
	{
		if (cell.startsWith("C"))
		{
			if (gene.charAt(1) == cell.charAt(1))
			{
				return generateValue(HIGH);
			}
			else
			{
				return generateValue(LOW);
			}
		}
		else if (gene.startsWith("G3"))
		{
			return generateValue(LOW, SD_BRIDGE);
		}
		else
		{
			int steps = C1_C2_BRIDGE_SIZE - 1;
			double pos = Double.parseDouble(cell.substring(2));
			if (gene.startsWith("G1")) pos = steps - pos;

			double target = LOW + ((HIGH - LOW) * (pos / steps));
			return generateValue(target, SD_BRIDGE);
		}
	}

	static double generateValue(double target)
	{
		return generateValue(target, SD);
	}

	static double generateValue(double target, double sd)
	{
		double value = target + (rand.nextGaussian() * sd);
		if (value < 0) value = 0;
		return value;
	}
}