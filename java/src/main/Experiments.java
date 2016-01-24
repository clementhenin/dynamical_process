package main;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

public class Experiments {

	/**
	 * First experiment, time evolution of concentration of C for several values of D
	 */
	public static void experiment1() {
		System.out.println("================Experiment 1=================") ;
		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*36 ;
		int nbB = 4500*36 ;
		int gridSize = 800*6 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double ro = 0 ;
		
		int nb_step = 10000 ;
		int time_step = 1 ;
		
		//Ranged parameters
		double[] D_range = new double[6] ;
		
		String[] labels = new String[D_range.length +1] ;
		labels[0] = "nb_step" ;
		
		for(int i = 0 ; i < D_range.length ; i ++)
		{
			double D = Math.pow(10, -i/2.) ;
			D_range[i] = D ;
			labels[i+1] = "D=" + D ;
		}

		//Store the resulting values
		double[][] values = new double[nb_step][D_range.length + 1] ;
		
		//Simulation loop
		int k = 1 ;
		for(double D0 : D_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pon = pon ;
			bs.poff = poff ;

			for(int i =0 ; i < nb_step ; i ++)
			{
				values[i][0] = bs.nb_simul * dt ;
				values[i][k] = bs.getC()/(double)nbA ;
				bs.simul(time_step) ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment1.csv", labels) ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}
	
	/**
	 * Time evolution of concentration of C for several values of ro
	 */
	public static void experiment2(){
		System.out.println("================Experiment 2=================") ;
 		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*36 ;
		int nbB = 4500*36 ;
		int gridSize = 800*6 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double D0 = 1 ;
		
		int nb_step = 10000 ;
		int time_step = 1 ;
		
		//Ranged parameters
		double[] ro_range = new double[8] ;
		
		String[] labels = new String[ro_range.length +1] ;
		labels[0] = "nb_step" ;
		
		for(int i = 0 ; i < ro_range.length ; i ++)
		{
			double ro = 0.05 * i ;
			ro_range[i] = ro ;
			labels[i+1] = "ro=" + ro ;
		}

		//Store the resulting values
		double[][] values = new double[nb_step][ro_range.length + 1] ;
		
		//Simulation loop
		int k = 1 ;
		for(double ro : ro_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pon = pon ;
			bs.poff = poff ;

			for(int i =0 ; i < nb_step ; i ++)
			{
				values[i][0] = bs.nb_simul * dt ;
				values[i][k] = bs.getC()/(double)nbA ;
				bs.simul(time_step) ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment2.csv", labels) ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
		
	}

	
	/**
	 * Equilibrium concentration of C as a function of D
	 */
	public static void experiment3() {
		System.out.println("================Experiment 3=================") ;
 		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*36 ;
		int nbB = 4500*36 ;
		int gridSize = 800*6 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double ro = 0 ;
		
		int nb_step = 10000 ;
		int avg_window = 200 ;
		
		//Ranged parameters
		double[] D_range = new double[12] ;
		
		for(int i = 0 ; i < D_range.length ; i ++)
		{
			double D = (i+1) / ((double) D_range.length) ;
			D_range[i] = D ;
		}

		//Store the resulting values
		double[][] values = new double[D_range.length][2] ;
		
		//Simulation loop
		int k = 0 ;
		for(double D0 : D_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pon = pon ;
			bs.poff = poff ;

			bs.simul(nb_step) ;
			
			values[k][0] = D0 ;
			for(int i =0 ; i < avg_window ; i ++)
			{
				values[k][1] += bs.getC()/(avg_window * (double)nbA) ;
				bs.simul() ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment3.csv", "D", "ratio C/C_max") ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}
	

	/**
	 * Equilibrium concentration of C as a function of ro
	 */
	public static void experiment4() {
		System.out.println("================Experiment 4=================") ;
 		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*16 ;
		int nbB = 4500*16 ;
		int gridSize = 800*4 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double D0 = 1 ;
		
		int nb_step = 10000 ;
		int avg_window = 100 ;
		
		//Ranged parameters
		double[] ro_range = new double[10] ;
		
		for(int i = 0 ; i < ro_range.length ; i ++)
		{
			double ro = (i+1) / ((double) ro_range.length)*0.4 ;
			ro_range[i] = ro ;
		}

		//Store the resulting values
		double[][] values = new double[ro_range.length][2] ;
		
		//Simulation loop
		int k = 0 ;
		for(double ro : ro_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pon = pon ;
			bs.poff = poff ;

			bs.simul(nb_step) ;
			
			values[k][0] = ro ;
			for(int i =0 ; i < avg_window ; i ++)
			{
				values[k][1] += bs.getC()/(avg_window * (double)nbA) ;
				bs.simul() ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment4.csv", "ro", "ratio C/C_max") ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}

	/**
	 * Equilibrium concentration of C as a function of D1
	 */
	public static void experiment5() {
		System.out.println("================Experiment 5=================") ;
 		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*9 ;
		int nbB = 4500*9 ;
		int gridSize = 800*3 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double D0 = 1 ;
		double ro = 0.35 ;
		
		int nb_step = 10000 ;
		int avg_window = 100 ;
		
		//Ranged parameters
		double[] D1_range = new double[20] ;
		
		for(int i = 0 ; i < D1_range.length ; i ++)
		{
			double D1 = (i+1) / ((double) D1_range.length)*D0 ;
			D1_range[i] = D1 ;
		}

		//Store the resulting values
		double[][] values = new double[D1_range.length][2] ;
		
		//Simulation loop
		int k = 0 ;
		for(double D1 : D1_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pmove1 = 4 * dt / (dx*dx) * D1 ;
			bs.patchSize = gridSize / 2 ;
			
			bs.pon = pon ;
			bs.poff = poff ;

			bs.simul(nb_step) ;
			
			values[k][0] = D1/D0 ;
			for(int i =0 ; i < avg_window ; i ++)
			{
				values[k][1] += bs.getC()/(avg_window * (double)nbA) ;
				bs.simul() ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment5.csv", "D1/D0", "ratio C/C_max") ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}

	/**
	 * Equilibrium concentration of C as a function of phi (fraction of patch area) 
	 */
	public static void experiment6() {
		System.out.println("================Experiment 6=================") ;
 		//Fixed parameters
		double dt = 1 ;
		double dx = 2 ; 
		int nbA = 100*9 ;
		int nbB = 4500*9 ;
		int gridSize = 800*3 ;
		double pon = 0.1 ; 
		double poff = 1e-3 ;
		double D0 = 1 ;
		double D1 = 0.5 ;
		double ro = 0.35 ;
		
		int nb_step = 10000 ;
		int avg_window = 100 ;
		
		//Ranged parameters
		double[] phi_range = new double[20] ;
		
		for(int i = 0 ; i < phi_range.length ; i ++)
		{
			double phi = (i+1) / ((double) phi_range.length) ;
			phi_range[i] = phi ;
		}

		//Store the resulting values
		double[][] values = new double[phi_range.length][2] ;
		
		//Simulation loop
		int k = 0 ;
		for(double phi : phi_range)
		{
			System.out.println("loop : " + k) ;
			BrownianSimulation bs = new BrownianSimulation(gridSize, nbA, nbB, ro) ;
			bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
			bs.pmove1 = 4 * dt / (dx*dx) * D1 ;
			bs.patchSize = (int) (gridSize * Math.sqrt(phi)) ;
			
			bs.pon = pon ;
			bs.poff = poff ;

			bs.simul(nb_step) ;
			
			values[k][0] = phi ;
			for(int i =0 ; i < avg_window ; i ++)
			{
				values[k][1] += bs.getC()/(avg_window * (double)nbA) ;
				bs.simul() ;
			}
			k ++ ;
		}

		try {
			//Write output
			BrownianSimulation.writeCsv(values, "experiment6.csv", "phi", "ratio C/C_max") ;
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}
	
	
	
	public static void main(String[] args) {
		//experiment1() ;
		experiment2() ;
		//experiment3() ;
		//experiment4() ;
		//experiment5() ;
		//experiment6() ;
		
		//TODO : refaire 4 avec grille plus grande
		//TODO : redo 5 and 6
	}

}
