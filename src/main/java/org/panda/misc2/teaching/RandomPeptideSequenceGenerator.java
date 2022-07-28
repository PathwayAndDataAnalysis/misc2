package org.panda.misc2.teaching;

import java.util.Random;

public class RandomPeptideSequenceGenerator
{
	public static final String AA = "GAVLIDENQPFWKCMYRHST";
	public static final Random RAND = new Random();

	public static void main(String[] args)
	{
		printAround("S", 5, 20);
	}

	public static void printAround(String centerAA, int halfL, int num)
	{
		for (int i = 0; i < num; i++)
		{
			System.out.println(getARandomSeq(halfL) + centerAA + getARandomSeq(halfL));
		}
	}

	private static String getARandomSeq(int length)
	{
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < length; i++)
		{
			sb.append(AA.charAt(RAND.nextInt(AA.length())));
		}
		return sb.toString();
	}
}
