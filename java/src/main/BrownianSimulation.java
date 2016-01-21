package main;

import java.awt.Point;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Class used for simulating Brownian motion of the three types of 
 * particles and their reactions 
 * @author marc
 *
 */
public class BrownianSimulation {

	Random rand = new Random() ;

	public double pon = 0.01 ; //Probability of A+B->C reacting 
	public double poff = 1e-3; // Probability of C->A+B decomposing

	
	public double pmove0 = 0.5 ; //Probability of moving inside the patch
	public double pmove1 = 0.5 ; //Probability of moving outside the patch
	
	public int patchSize ; //Size of the patch

	public int nb_simul = 0 ; // Number of simulation steps up to now

	
	private int nb_C = 0 ; //Number of C particles in the system
	List<Point> positions = new ArrayList<>() ; //List of position where the particles are

	private final int gridSize ; //Size of the grid for the simulation
	int[][] aParticles ; //Table containing the number of A particles at position (i,j)
	int[][] bParticles ; //Same with B particles
	int[][] cParticles ; //Same with C particles
	
	//Additional tables used for the simulations
	private int[][] aPartNew ; 
	private int[][] bPartNew ;
	private int[][] cPartNew ;
	final boolean[][] obstacles ; // table containing the obstacles
	final int[][] visited ; //table giving the last time step a given cell was visited to remove duplicates from 'positions'

	public BrownianSimulation(int gridSize, int nb_A, int nb_B, double ro) {
		
		//Initialize memory
		this.gridSize = gridSize ;
		aParticles = new int[gridSize][gridSize] ;
		bParticles = new int[gridSize][gridSize] ;
		cParticles = new int[gridSize][gridSize] ;
		aPartNew = new int[gridSize][gridSize] ;
		bPartNew = new int[gridSize][gridSize] ;
		cPartNew = new int[gridSize][gridSize] ;
		obstacles = new boolean[gridSize][gridSize] ;
		visited = new int[gridSize][gridSize] ;

		patchSize = gridSize ; 
		init(nb_A, nb_B, ro) ;
	}

	/**
	 * Initialize the position of the particles and obstacles
	 * @param nbA
	 * @param nbB
	 * @param ro
	 */
	private void init(int nbA, int nbB, double ro) {
		for(int i = 0 ; i < gridSize ; i++)
		{
			for(int j = 0 ; j < gridSize ; j++)
			{
				if(rand.nextFloat() < ro)
					obstacles[i][j] = true ;
				else 
					obstacles[i][j] = false ;
			}
		}

		int i = 0 ;
		while(i < nbA)
		{
			int x = rand.nextInt(gridSize) ;
			int y = rand.nextInt(gridSize) ;
			if(!obstacles[x][y])
			{
				i+= 1 ;
				aParticles[x][y] += 1 ;
				if(visited[x][y] != -1)
					positions.add(new Point(x,y)) ;
				visited[x][y] = -1 ;
			}
		}
		i = 0 ;
		while(i < nbB)
		{
			int x = rand.nextInt(gridSize) ;
			int y = rand.nextInt(gridSize) ;
			if(!obstacles[x][y])
			{
				i+= 1 ;
				bParticles[x][y] += 1 ;
				if(visited[x][y] != -1)
					positions.add(new Point(x,y)) ;
				visited[x][y] = -1 ;
			}
		}
	}

	
	/**
	 * Returns the probability of a particle moving depending of the position
	 * @param p : position (i,j)
	 * @return pmove(i,j) = pmove0 if in the patch, pmove1 otherwise
	 */
	private double getPmove(Point p){
		int border = (gridSize - patchSize) /2 ; 
		if(p.x > border && p.x < gridSize - border && p.y > border && p.y < gridSize - border)
			return pmove0;
		else 
			return pmove1 ;
	}
	
	/**
	 * @return Number of C particles in the system
	 */
	public int getC() {
		return nb_C ;
	}

	/**
	 * Runs one simulation step
	 */
	public void simul(){
		nb_simul += 1 ;
		List<Point> newPositions = new ArrayList<>() ;

		for(Point p : positions)
		{
			reactAB(p) ;
			reactC(p) ;
		}

		for(Point p : positions)
		{
			move(aParticles, aPartNew, newPositions, p) ;
			move(bParticles, bPartNew, newPositions, p) ;
			move(cParticles, cPartNew, newPositions, p) ;
		}
		swap() ;
		positions = newPositions ;
	}

	/**
	 * Swaps the tables aParticles/aNewPart, bPartices/....
	 */
	private void swap(){
		int[][] temp ;
		temp = aPartNew ;
		aPartNew = aParticles ;
		aParticles = temp ;

		temp = bPartNew ;
		bPartNew = bParticles ;
		bParticles = temp ;

		temp = cPartNew ;
		cPartNew = cParticles ;
		cParticles = temp ;		
	}

	/**
	 * Simulates for the given number of steps
	 * @param nbStep : number of steps for the simulation
	 */
	public void simul(int nbStep) {
		for(int i = 0 ; i < nbStep ; i++)
		{
			simul() ;
		}
	}

	/**
	 * Simulates for the given amount of time
	 * @param sec : time in second
	 */
	public void simul(double sec){
		double start = System.currentTimeMillis() ;
		while(System.currentTimeMillis() - start < sec * 1000)
		{
			simul() ;
		}
	}

	/**
	 * goes through the table and performs for C->A+B reactions
	 * @param p
	 */
	private void reactC(Point p) {
		int nb_react = 0 ;
		for(int i = 0 ; i < cParticles[p.x][p.y] ; i++)
		{
			if(rand.nextFloat() < poff)
				nb_react += 1;
		}
		cParticles[p.x][p.y] -= nb_react ;
		aParticles[p.x][p.y] += nb_react ;
		bParticles[p.x][p.y] += nb_react ;
		nb_C -= nb_react ;

	}

	/**
	 * goes through the table and performs for A+B->C reactions
	 * @param p
	 */
	private void reactAB(Point p) {
		if(aParticles[p.x][p.y] == 0 || bParticles[p.x][p.y] == 0)
			return ;
		int nb_react = 0 ;
		for(int i = 0 ; i < aParticles[p.x][p.y]; i++)
		{
			for(int j = 0 ; j < bParticles[p.x][p.y] - nb_react ; j++)
			{

				if(rand.nextFloat() < pon)
				{
					nb_react += 1 ;
					break ;
				}
			}
		}
		cParticles[p.x][p.y] += nb_react ;
		aParticles[p.x][p.y] -= nb_react ;
		bParticles[p.x][p.y] -= nb_react ;
		nb_C += nb_react ;
	}

	/**
	 * Fits the given position inside the grid
	 * @param x
	 * @param y
	 * @return
	 */
	private Point fitGrid(int x, int y){
		return new Point(
				Math.max(0, Math.min(x, gridSize -1)), 
				Math.max(0, Math.min(y, gridSize -1))) ;
	}

	/**
	 * Generic function for moving particles for one time step
	 * @param source
	 * @param target
	 * @param newPositions
	 * @param p
	 */
	private void move(int[][] source, int[][] target, List<Point> newPositions, Point p){
		int[] newVals = new int[5] ;

		double pmove = getPmove(p) ;
		
		for(int i = 0 ; i < source[p.x][p.y] ; i++)
		{
			if(rand.nextFloat() < pmove)
				newVals[rand.nextInt(4) + 1] += 1 ;
			else 
				newVals[0] += 1 ;
		}
		if(newVals[0] != 0)
		{
			target[p.x][p.y] += newVals[0] ;
			if(visited[p.x][p.y] != nb_simul)
			{
				visited[p.x][p.y] = nb_simul ;
				newPositions.add(p) ;
			}
		}
		if(newVals[1] != 0)
		{
			Point p2 = fitGrid(p.x + 1, p.y) ;
			target[p2.x][p2.y] += newVals[1] ;
			if(visited[p2.x][p2.y] != nb_simul)
			{
				visited[p2.x][p2.y] = nb_simul ;
				newPositions.add(p2) ;
			}
		}
		if(newVals[2] != 0)
		{
			Point p2 = fitGrid(p.x - 1, p.y) ;
			target[p2.x][p2.y] += newVals[2] ;
			if(visited[p2.x][p2.y] != nb_simul)
			{
				visited[p2.x][p2.y] = nb_simul ;
				newPositions.add(p2) ;
			}
		}
		if(newVals[3] != 0)
		{
			Point p2 = fitGrid(p.x, p.y+1) ;
			target[p2.x][p2.y] += newVals[3] ;
			if(visited[p2.x][p2.y] != nb_simul)
			{
				visited[p2.x][p2.y] = nb_simul ;
				newPositions.add(p2) ;
			}
		}
		if(newVals[4] != 0)
		{
			Point p2 = fitGrid(p.x, p.y-1) ;
			target[p2.x][p2.y] += newVals[4] ;
			if(visited[p2.x][p2.y] != nb_simul)
			{
				visited[p2.x][p2.y] = nb_simul ;
				newPositions.add(p2) ;
			}
		}
		source[p.x][p.y] = 0 ;
	}

	/**
	 * Utility function for writing data to a CSV file
	 * @param values : the data
	 * @param fileName : the name of the file
	 * @param labels : the name of each column
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 */
	public static void writeCsv(double[][] values, String fileName, String... labels) throws FileNotFoundException, UnsupportedEncodingException {
		if(values[0].length != labels.length)
		{
			System.out.println("Incorrect number of labels for the columns name") ;
			return ;
		}
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");
		for(int i = 0 ; i < labels.length - 1 ; i++)
			writer.print(labels[i] + "\t") ;
		writer.println(labels[labels.length - 1]) ;
		for(int i = 0 ; i < values.length ; i++)
		{
			for(int j = 0 ; j < values[i].length -1 ; j++)
				writer.print(values[i][j] + "\t") ;
			writer.println(values[i][values[i].length -1]) ;
		}
		writer.close();
	}

	/**
	 * Performs check on the state of the simulation to verify that there were no errors during the simulation
	 */
	public void test(){
		for(int i = 0 ; i < gridSize ; i++)
			for(int j = 0 ; j< gridSize ; j++)
			{
				if(aParticles[i][j] < 0)
					System.out.println("Negative value for A: " + i + "," + j) ;

				if(bParticles[i][j] < 0)
					System.out.println("Negative value for B: " + i + "," + j) ;


				if(cParticles[i][j] < 0)
					System.out.println("Negative value for C: " + i + "," + j) ;

				if(aParticles[i][j] != 0 || bParticles[i][j] != 0 || cParticles[i][j] != 0)
					if(visited[i][j] != nb_simul)
					{
						System.out.println("Lost position : " + i + "," + j) ;
						System.out.println("Particles : " + aParticles[i][j] + "," + bParticles[i][j]  +"," + cParticles[i][j]) ;
					}	

				if(visited[i][j] == nb_simul)
				{
					if(aParticles[i][j] == 0 && bParticles[i][j] == 0 && cParticles[i][j] == 0)
						System.out.println("Added position : " + i + "," + j) ;
				}

			}
	}

	public static void main(String[] args) {
		double dt = 1 ; 
		double dx = 2 ;
		double D0 = 1 ;
		int nbA = 100 ;
		BrownianSimulation bs = new BrownianSimulation(800, nbA, 4500, 0.0) ;
		bs.pmove0 = 4 * dt / (dx*dx) * D0 ;
		bs.pon = 0.1 ;
		bs.poff = 1e-3 ;


		int nb_step = 30000 ;
		int step = 1 ;
		double[][] values = new double[nb_step][2] ;
		for(int i =0 ; i < nb_step ; i ++)
		{
			values[i][0] = bs.nb_simul ;
			values[i][1] = bs.getC()/(double)nbA ;
			bs.simul(step) ;
		}
		try {
			writeCsv(values, "output.csv", "nb_step", "C") ;


		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		System.out.println("............................DONE") ;
	}

}
