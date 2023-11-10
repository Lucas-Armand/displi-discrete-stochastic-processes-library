# DiSPLi: Discrete Stochastic Processes Library for Python

Displi is a Python library designed for discrete stochastic processes, providing analytical solutions for queueing systems based on queueing theory. It includes implementations for M/M/S/K and M/G/1/K queueing models, allowing users to analyze and understand the behavior of queuing systems.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Installation

To install Displi, you can use `pip`:

```bash
pip install displi
```
## Usage

Here's a quick overview of the available queueing models:
MMSK (M/M/S/K)

```python
from displi import MMSK

# Create an MMSK queueing system
mmsk = MMSK(arrival_rate, service_rate, servers, capacity)

# Calculate traffic intensity
ro = mmsk.calculate_traffic_intensity()

# Calculate probability of having zero customers
p0 = mmsk.calculate_probability_zero_customers()

# ... (and more)
```

MG1K (M/G/1/K)

```python
from displi import MG1K

# Create an MG1K queueing system
mg1k = MG1K(arrival_rate, service_rate, variance, capacity)

# Calculate traffic intensity
ro = mg1k.calculate_traffic_intensity()

# Calculate probability of having K or more customers
pk = mg1k.calculate_probability_k()

# ... (and more)
```


## Examples

Here's a simple example of using the MMSK model:

```python

from displi import MMSK

# Create an MMSK queueing system
mmsk = MMSK(arrival_rate=5, service_rate=8, servers=3, capacity=10)

# Calculate traffic intensity
ro = mmsk.calculate_traffic_intensity()

# Print the result
print(f"Traffic Intensity (ro): {ro}")
```

## Documentation

For detailed documentation and additional details about the methods and parameters, refer to the Documentation section.
Contributing

If you want to contribute to the development of Displi, feel free to submit issues or pull requests. Your contributions are welcome!
## License

Displi is distributed under the MIT License.
