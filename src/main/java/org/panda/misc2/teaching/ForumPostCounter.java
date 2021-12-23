package org.panda.misc2.teaching;

import org.panda.utility.TermCounter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * On the discussion board, select all, press collect, press Print Preview, on the popup cancel the print dialog, save
 * as single file.
 */
public class ForumPostCounter
{
	public static void main(String[] args) throws FileNotFoundException
	{
		count("/home/ozgunbabur/Downloads/temp.html");
	}

	public static void count(String file) throws FileNotFoundException
	{
		TermCounter tc = new TermCounter();
		Scanner sc = new Scanner(new File(file));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			if (line.equals("     =20"))
			{
				String nameline = sc.nextLine();
				String underLine = sc.nextLine();

				if (underLine.equals(" =20"))
				{
					tc.addTerm(nameline.trim());
				}
			}
		}
		sc.close();

		tc.print();
	}
}
