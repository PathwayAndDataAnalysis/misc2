package org.panda.misc2.teaching;

import org.panda.utility.ArrayUtil;

import java.util.Arrays;
import java.util.Random;

public class DataGenerator
{
	public static void main(String[] args)
	{
//		printNoisyLine(1, 2, 0.6, 100, -5, 5);
//		printNoisyLine(new double[]{5, 23, -6, 8}, new double[]{1, 30, 10, 50}, new double[]{1, 100, -20, 50}, 100, 0.05);
//		printNoisyLine(new double[]{4, 2, -1, 0, -3}, new double[]{1, 3, 5, 0.5, 2}, new double[]{10, 5, -10, -5, 0}, 200, 0.05);
//		printLogisticReg2DDatasetwithMultFeat(new double[]{-15, 0, 0, 1, }, new double[]{0, 10}, new double[]{0, 10}, 100);
//		printLogisticRegDataset(new double[]{3, 1, -3, 2}, new double[][]{new double[]{0, 10}, new double[]{5, 10}, new double[]{-3, 15}},  50);
//		printLogisticReg2DDatasetWithAnchor(new double[]{-2, 2}, new double[]{-2, 2}, 50, new double[][]{new double[]{0, 0}}, new double[]{1.6});
		print2DLogisticRegDataset3Classes(new double[]{-2, 2}, new double[]{-2, 2}, 100, new double[][]{new double[]{-1.5, 0}, new double[]{-1, 1}, new double[]{-1, -1}});
	}

	public static void printNoisyLine(double t0, double t1, double sd, int m, double xMin, double xMax)
	{
		Random rand = new Random();
		double range = xMax - xMin;

		for (int i = 0; i < m; i++)
		{
			double x = (rand.nextDouble() * range) + xMin;
			double y = t0 + (t1 * x);

			x += rand.nextGaussian()*sd;
			y += rand.nextGaussian()*sd;

			System.out.println(x + "," + y);
		}
	}

	public static void printNoisyLine(double[] theta, double[] sd, double[] mean, int m, double noiseSDFactor)
	{
		if (theta.length != sd.length || sd.length != mean.length)
		{
			throw new IllegalArgumentException("Wrong sizes!");
		}

		Random rand = new Random();

		for (int i = 0; i < m; i++)
		{
			double[] x = new double[theta.length];
			x[0] = 1;
			for (int j = 1; j < x.length; j++)
			{
				x[j] = (rand.nextGaussian() * sd[j]) + mean[j];
			}
			double y = innerProduct(theta, x);

			for (int j = 1; j < x.length; j++)
			{
				x[j] += rand.nextGaussian()*(sd[j] * noiseSDFactor);
			}
			y += rand.nextGaussian()*noiseSDFactor;

			for (int j = 1; j < x.length; j++)
			{
				System.out.print(x[j] + ",");
			}
			System.out.println(y);
		}

	}

	public static double innerProduct(double[] v1, double[] v2)
	{
		double sum = 0;
		for (int i = 0; i < v1.length; i++)
		{
			sum += v1[i] * v2[i];
		}
		return sum;
	}

	public static void printLogisticReg2DDatasetwithMultFeat(double[] theta, double[] range1, double[] range2, int m)
	{
		Random  rand = new Random();
		int[] cnt = new int[]{0, 0};
		for (int i = 0; i < m; i++)
		{
			double[] x = new double[theta.length];
			x[0] = 1;

			x[1] = randomInRange(rand, range1);
			x[2] = randomInRange(rand, range2);
			if (theta.length > 3) x[3] = x[1] * x[2];

			int y = innerProduct(theta, x) < 0 ? 0 : 1;
			cnt[y]++;

			System.out.println(x[1] + "," + x[2] + "," + y);
		}

		System.out.println("cnt = " + Arrays.toString(cnt));
	}

	public static void printLogisticReg2DDatasetWithAnchor(double[] range1, double[] range2, int m, double[][] anchors, double[] proximity)
	{
		Random  rand = new Random();
		int[] cnt = new int[]{0, 0};
		for (int i = 0; i < m; i++)
		{
			double[] x = new double[3];
			x[0] = 1;

			x[1] = randomInRange(rand, range1);
			x[2] = randomInRange(rand, range2);

			boolean qualifies = false;
			for (int j = 0; j < anchors.length; j++)
			{
				double[] anchor = anchors[j];
				double dist = Math.sqrt(Math.pow(x[1] - anchor[0], 2) + Math.pow(x[2] - anchor[1], 2));
				qualifies = dist < proximity[j];
				if (qualifies) break;
			}

			int y = qualifies && Math.random() < 0.7 ? 1 : 0;
			cnt[y]++;

			System.out.println(x[1] + "," + x[2] + "," + y);
		}

		System.out.println("cnt = " + Arrays.toString(cnt));
	}

	public static void printLogisticRegDataset(double[] theta, double[][] ranges, int m)
	{
		Random  rand = new Random();
		int[] cnt = new int[]{0, 0};
		for (int i = 0; i < m; i++)
		{
			double[] x = new double[theta.length];
			x[0] = 1;

			for (int j = 1; j < theta.length; j++)
			{
				x[j] = randomInRange(rand, ranges[j-1]);
			}

			int y = innerProduct(theta, x) < 0 ? 0 : 1;
			cnt[y]++;

			for (int j = 1; j < x.length; j++)
			{
				System.out.print(x[j] + ",");
			}
			System.out.println(y);
		}

		System.out.println("cnt = " + Arrays.toString(cnt));
	}

	public static void print2DLogisticRegDataset3Classes(double[] range1, double[] range2, int m, double[][] anchors)
	{
		Random  rand = new Random();
		int[] cnt = new int[anchors.length];

		for (int i = 0; i < m; i++)
		{
			double[] x = new double[3];
			x[0] = 1;

			x[1] = randomInRange(rand, range1);
			x[2] = randomInRange(rand, range2);

			double[] dist = new double[anchors.length];

			for (int j = 0; j < anchors.length; j++)
			{
				double[] anchor = anchors[j];
				dist[j] = Math.sqrt(Math.pow(x[1] - anchor[0], 2) + Math.pow(x[2] - anchor[1], 2));
			}

			int y = ArrayUtil.indexOfMax(dist) + 1;
			cnt[y-1]++;

			System.out.println(x[1] + "," + x[2] + "," + y);
		}

		System.out.println("cnt = " + Arrays.toString(cnt));
	}

	public static double randomInRange(Random rand, double[] range)
	{
		double r = rand.nextDouble();
		double width = range[1] - range[0];
		return (r * width) + range[0];
	}
}
