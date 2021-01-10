#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "neuron.h"

neuron *create_neuron(int num_inputs, int num_outputs){
	neuron *output;

	output = malloc(sizeof(neuron));
	output->num_inputs = num_inputs;
	output->num_outputs = num_outputs;
	output->input_neurons = malloc(sizeof(neuron *)*num_inputs);
	output->output_neurons = malloc(sizeof(neuron *)*num_outputs);
	output->weights = calloc(num_inputs, sizeof(double));
	output->output = 0;
	output->next_output = 0;

	return output;
}

static void resize_outputs(neuron *n){
	if(!(n->output_neurons = realloc(n->output_neurons, sizeof(neuron *)*n->num_outputs))){
		fprintf(stderr, "Error: out of memory.");
		exit(1);
	}
}

static void resize_inputs(neuron *n){
	if(!(n->input_neurons = realloc(n->input_neurons, sizeof(neuron *)*n->num_inputs))){
		fprintf(stderr, "Error: out of memory.");
		exit(1);
	}
	if(!(n->weights = realloc(n->weights, sizeof(double)*n->num_inputs))){
		fprintf(stderr, "Error: out of memory.");
		exit(1);
	}
}

void connect_neurons(neuron *from, neuron *to, double weight){
	from->num_outputs++;
	resize_outputs(from);
	from->output_neurons[from->num_outputs - 1] = to;

	to->num_inputs++;
	resize_inputs(to);
	to->input_neurons[to->num_inputs - 1] = from;
	to->weights[to->num_inputs - 1] = weight;
}

double update_neuron(neuron *n){
	int i;
	double sum = 0;

	for(i = 0; i < n->num_inputs; i++){
		sum += n->input_neurons[i]->output*n->weights[i];
	}

	n->next_output = 1/(1 + exp(-sum));

	return (n->output - n->next_output)*(n->output - n->next_output);
}

void update_neurons(neuron *neurons, int num_neurons, int count, double min_change){
	int i;
	int j;
	double difference;
	double biggest_diff;

	for(i = 0; i < count; i++){
		biggest_diff = -1;
		for(j = 0; j < num_neurons; j++){
			difference = update_neuron(neurons + j);
			if(difference > biggest_diff){
				biggest_diff = difference;
			}
		}
		for(j = 0; j < num_neurons; j++){
			neurons[j].output = neurons[j].next_output;
		}
		if(biggest_diff < min_change){
			return;
		}
	}
}

