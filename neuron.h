typedef struct neuron neuron;

struct neuron{
	int num_inputs;
	int num_outputs;
	double output;
	double next_output;
	neuron **input_neurons;
	neuron **output_neurons;
	double *weights;
};

neuron *create_neuron(int num_inputs, int num_outputs);

void connect_neurons(neuron *from, neuron *to, double weight);

double update_neuron(neuron *n);

void update_neurons(neuron *neurons, int num_neurons, int count, double min_change);
