#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "neuron.h"
#include "neat.h"

double add_connection_chance = 0.3;
double mutate_chance = 0.05;
double split_chance = 0.05;
double toggle_activation_chance = 0.05;
double random_weight_magnitude = 2;

double disjoint_distance_coef = 0.1;
double weight_diff_coef = 0.5;

long int global_innovation_num = 0;

int num_species;
species **current_species;

static double random_double(){
	return ((double) rand())/RAND_MAX;
}

static double get_random_weight(){
	double output;
	output = -log(1/random_double() - 1)*random_weight_magnitude;
	if(isnan(output)){
		return 1;
	}
	return output;
}

genome *create_genome(int num_neurons, int num_genes, int species){
	genome *output;

	output = malloc(sizeof(genome));
	output->num_genes = num_genes;
	output->num_neurons = num_neurons;
	output->genes = malloc(sizeof(gene)*num_genes);
	output->species = species;

	return output;
}

gene create_gene(int from, int to, double weight, int innovation_num){
	return (gene) {.from = from, .to = to, .weight = weight, .innovation_num = innovation_num, .activated = 1};
}

void add_gene(genome *g, gene new){
	g->num_genes++;
	if(!(g->genes = realloc(g->genes, sizeof(gene)*g->num_genes))){
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
	g->genes[g->num_genes - 1] = new;
}

void mutate(genome *g){
	int i;
	int orig_num_genes;

	orig_num_genes = g->num_genes;

	for(i = 0; i < orig_num_genes; i++){
		if(random_double() < toggle_activation_chance){
			g->genes[i].activated = !g->genes[i].activated;
		}
		if(random_double() < mutate_chance){
			g->genes[i].weight = get_random_weight();
		}
		if(random_double() < split_chance){
			g->genes[i].activated = 0;
			g->num_neurons++;
			add_gene(g, (gene) {.from = g->genes[i].from, .to = g->num_neurons - 1, .weight = 1, .innovation_num = global_innovation_num, .activated = 1});
			global_innovation_num++;
			add_gene(g, (gene) {.from = g->num_neurons - 1, .to = g->genes[i].to, .weight = g->genes[i].weight, .innovation_num = global_innovation_num, .activated = 1});
			global_innovation_num++;
		}
	}

	if(random_double() < add_connection_chance){
		add_gene(g, (gene) {.from = rand()%g->num_neurons, .to = rand()%g->num_neurons, .weight = get_random_weight(), .innovation_num = global_innovation_num, .activated = 1});
		global_innovation_num++;
	}
}

genome *mate(genome *more_fit, genome *less_fit){
	genome *output;
	int gene_index0;
	int gene_index1 = 0;

	output = create_genome(more_fit->num_neurons, more_fit->num_genes, more_fit->species);

	for(gene_index0 = 0; gene_index0 < more_fit->num_genes; gene_index0++){
		while(gene_index1 < less_fit->num_genes && more_fit->genes[gene_index0].innovation_num > less_fit->genes[gene_index1].innovation_num){
			gene_index1++;
		}
		if(gene_index1 < less_fit->num_genes && less_fit->genes[gene_index1].innovation_num == more_fit->genes[gene_index0].innovation_num && rand()%2){
			output->genes[gene_index0] = less_fit->genes[gene_index1];
		} else {
			output->genes[gene_index0] = more_fit->genes[gene_index0];
		}
	}

	return output;
}

double genome_distance(genome *a, genome *b){
	double total_weight_distance;
	int num_shared = 0;
	int num_disjoint = 0;
	int gene_index0;
	int gene_index1;

	for(gene_index0 = 0; gene_index0 < a->num_genes; gene_index0++){
		while(gene_index1 < b->num_genes && a->genes[gene_index0].innovation_num > b->genes[gene_index1].innovation_num){
			gene_index1++;
			num_disjoint++;
		}
		if(gene_index1 < b->num_genes && a->genes[gene_index0].innovation_num == b->genes[gene_index1].innovation_num){
			gene_index1++;
			num_shared++;
			total_weight_distance += fabs(a->genes[gene_index0].weight - b->genes[gene_index1].weight);
		} else {
			num_disjoint++;
		}
	}

	return num_disjoint*disjoint_distance_coef + total_weight_distance/num_shared*weight_diff_coef;
}

organism *create_organism_from_genome(genome *g){
	organism *output;
	int i;

	output = malloc(sizeof(organism));
	output->num_neurons = g->num_neurons;
	output->neurons = malloc(sizeof(neuron)*g->num_neurons);
	for(i = 0; i < g->num_neurons; i++){
		output->neurons[i] = create_neuron(0, 0);
	}
	for(i = 0; i < g->num_genes; i++){
		if(g->genes[i].activated){
			connect_neurons(output->neurons[g->genes[i].from], output->neurons[g->genes[i].to], g->genes[i].weight);
		}
	}

	return output;
}

species *create_species(int num_genomes){
	species *output;

	output = malloc(sizeof(species));
	output->num_genomes = num_genomes;
	output->genomes = malloc(sizeof(genome *)*num_genomes);

	return output;
}

void add_genome(species *s, genome *new_genome){
	s->num_genomes++;
	if(!(s->genomes = realloc(s->genomes, sizeof(genome *)*s->num_genomes))){
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
	s->genomes[s->num_genomes - 1] = new_genome;
}

static void compare_genome(void *a, void *b){
	genome *genome_a;
	genome *genome_b;

	genome_a = a;
	genome_b = b;

	return genome_b->fitness - genome_a->fitness;
}

void sort_into_species(species ***s, int *num_species, genome **genomes, int num_genomes){
	int i;
	int j;

	*s = NULL;
	*num_species = 0;

	for(i = 0; i < num_genomes; i++){
		for(j = 0; j < *num_species; j++){
			if(genome_distance(genomes[i], (*s)[j]->genomes[0]) < distance_threshold){
				add_genome((*s)[j], genomes[i]);
				break;
			}
		}
		if(j == *num_species){
			++*num_species;
			*s = realloc(*s, sizeof(species *)*(*num_species));
			(*s)[*num_species - 1] = create_species(1);
			(*s)[*num_species - 1]->genomes[0] = genomes[i];
		}
	}
}

void neat(int num_organisms, int num_generations, int initial_neurons, double (*evaluate)(organism *, int, int)){
	genome **genomes;
	genome **next_genomes;
	species **species_list;
	organism *new_organism;
	int i;
	int j;
	int num_species;

	genomes = malloc(sizeof(genome *)*num_organisms);
	next_genomes = malloc(sizeof(genome *)*num_organisms);
	for(i = 0; i < num_organisms; i++){
		genomes[i] = create_genome(initial_neurons, 0, 0);
	}
	for(i = 0; i < num_generations; i++){
		for(j = 0; j < num_organisms; j++){
			mutate(genomes[j]);
			new_organism = create_organism_from_genome(genomes[j]);
			genomes[j]->fitness = evaluate(new_organism, i, genomes[j]->species);
		}
		sort_into_species(&species_list, &num_species, genomes, num_organisms);
		for(j = 0; j < num_species; j++){
			
		}
	}
	free(genomes);
	free(next_genomes);
}

int main(int argc, char **argv){
	genome *g;
	organism *o;

	g = create_genome(5, 0);
	add_gene(g, create_gene(0, 1, 1, 221));
	o = create_organism_from_genome(g);
	printf("hello world!\n");
}

