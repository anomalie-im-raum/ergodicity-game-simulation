import numpy as np
import matplotlib.pyplot as plt

class Ergodicity:

    def __init__(self, winning_factor, losing_factor):
        self.wealth_factors = [losing_factor, winning_factor]
        self.rng = np.random.default_rng()

    def tossing_coin(self):
        """Simulates the tossing of a fair coin.
    
        Returns:
            Either 0 for head or 1 for tail.
        """
        return self.rng.integers(0, 1, endpoint=True)

    def calculate_new_score(self, score):
        """Calculates the new score based on the previous score and the defined probabilities for head or tail.

        Args:
            score: The previous score.
    
        Returns:
            The new score determined by chance.
        """
        coin_toss = self.tossing_coin()
        #print("Coin toss", coin_toss)
        if coin_toss == 0:
            score *= self.wealth_factors[0]
        else:
            score *= self.wealth_factors[1]    
        return score

    def create_coin_toss_series(self, nr_coin_tosses):
        """Creates a series of score values.

        Args:
            score: The number of coin tosses determining the length of the series.
    
        Returns:
            A series of score values.
        """
        scores = []
        for i in range(nr_coin_tosses):
            #The initial score is set to 1
            previous_score = scores[i-1] if i > 0 else 1
            scores.append(self.calculate_new_score(previous_score)) 
        return scores

    def plot_score(self, time, scores, average_scores):
        """Creates two plots. The first depicting the scores for a single coin toss series.
           The second depicting the averaged score of the Monte-Carle simulation, simulating 
           many single coin tosses in parallel. 

        Args:
            time: The passed time, eg. the passed 52 weeks.
            scores: The scores for a single coin toss series. Single event.
            averaged_scores: The averaged socre of many coin toss series. Parallel events.
        """
        fig, ax = plt.subplots(1, 2)
        ax[1].plot(time, scores, label="Coin Toss")
        ax[1].set_title("Individual Score")
        ax[1].set_xlabel("Number of Coin Tosses")
        ax[1].set_ylabel("Score Multiple")
        ax[0].plot(time, average_scores, label="Trial Runs")
        ax[0].set_title("Ensemble Average Score")
        ax[0].set_xlabel("Number of Coin Tosses")
        ax[0].set_ylabel("Score Multiple")
        plt.show()

    def monte_carlo_simulation(self, nr_of_trials, nr_coin_tosses):
        """Creates the ensemble scores, consisting of the of Monte-Carle simulation results.

        Args:
            nr_of_trials: The number of times the simulation is run.
            nr_coin_tosses: The number of coin tosses for a single pass.
    
        Returns:
            The ensemble scores. 
        """
        ensemble_scores = []
        for i in range(nr_of_trials):
            scores = self.create_coin_toss_series(nr_coin_tosses)
            ensemble_scores.append(scores)
        return ensemble_scores

    def calculate_finite_ensemble_average(self, ensemble_scores, nr_coin_tosses, nr_of_trials):
        """Calculates the average over all simulations at a given time. The ensemble average is the average at a fixed time over many systems.
        Eg. the average score for week 1 over 1000 Monte-Carlo simulations.

        Args:
            ensemble_scores: The Monte-Carlo simulation results.
            nr_of_trials: The number of times the simulation is run.
            nr_coin_tosses: The number of coin tosses for a single pass.
    
        Returns:
            The calculated averages for all Monte-Carlo simulation for a given time. 
        """       
        average_score = []
        for i in range(nr_coin_tosses):
            sum = 0
            for (j, k) in enumerate(ensemble_scores):
                sum += ensemble_scores[j][i] 
                #print("Ensemble scores", j, k)
            average = sum / nr_of_trials
            average_score.append(average)
        return average_score

    def run_process_simulation(self, nr_coin_tosses, nr_of_trials):
        score = ergodicity.create_coin_toss_series(nr_coin_tosses)
        ensemble_score = ergodicity.monte_carlo_simulation(nr_of_trials, nr_coin_tosses)
        average_score = ergodicity.calculate_finite_ensemble_average(ensemble_score, nr_coin_tosses, nr_of_trials)
        ergodicity.plot_score(range(nr_coin_tosses), score, average_score)
     
ergodicity = Ergodicity(winning_factor=1.5, losing_factor=0.6)
ergodicity.run_process_simulation(nr_coin_tosses=52, nr_of_trials=1000)

     



