import math


class MMS:
    def __init__(self, arrival_rate, service_rate, servers):
        """
        Initializes an M/M/S queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of each server.
        :param servers: Number of servers.
        """
        self.arrival_rate = float(arrival_rate)
        self.service_rate = float(service_rate)
        self.servers = int(servers)
        self.name = "MMS"

    def calculate_traffic_intensity(self):
        """Calculates the traffic intensity (ro) of the queueing system."""
        return self.arrival_rate / (self.service_rate * self.servers)

    def calculate_probability_zero_customers(self):
        """Calculates the probability of having zero customers in the system (P0)."""
        arrival_rate, service_rate, servers = self.arrival_rate, self.service_rate, self.servers
        S1 = sum((arrival_rate / (service_rate * servers)) ** i / math.factorial(i) for i in range(servers))
        S2 = ((arrival_rate / (service_rate * servers)) ** servers) / (
                math.factorial(servers) * (1 - arrival_rate / (service_rate * servers)))
        return 1 / (S1 + S2)

    def calculate_probability_n_customers(self, n):
        """
        Calculates the probability of having 'n' customers in the system.

        :param n: Number of customers.
        """
        arrival_rate, service_rate, servers = self.arrival_rate, self.service_rate, self.servers
        p0 = self.calculate_probability_zero_customers()
        if n < servers:
            return ((arrival_rate / service_rate) ** n) * p0 / math.factorial(n)
        else:
            return ((arrival_rate / service_rate) ** n) * p0 / (math.factorial(servers) * (servers ** (n - servers)))

    def calculate_queue_probability(self):
        """Calculates the probability of a customer being in the queue (P_queue)."""
        s = self.servers
        S = sum(self.calculate_probability_n_customers(i) for i in range(s))
        return 1 - S

    def calculate_average_queue_length(self):
        """Calculates the average number of customers in the queue (Lq)."""
        arrival_rate, service_rate, servers, p0 = (
            self.arrival_rate, self.service_rate, self.servers, self.calculate_probability_zero_customers()
        )
        Lq = ((arrival_rate / service_rate) ** (servers + 1)) * p0 / (
                math.factorial(servers) * (1 - arrival_rate / (service_rate * servers)) ** 2)
        return Lq

    def calculate_average_system_length(self):
        """Calculates the average number of customers in the system (L)."""
        arrival_rate, service_rate, Lq = self.arrival_rate, self.service_rate, self.calculate_average_queue_length()
        L = Lq + arrival_rate / service_rate
        return L

    def calculate_system_time(self):
        """Calculates the average time a customer spends in the system (W)."""
        arrival_rate, L = self.arrival_rate, self.calculate_average_system_length()
        W = L / arrival_rate
        return W

    def calculate_queue_time(self):
        """Calculates the average time a customer spends in the queue (Wq)."""
        arrival_rate, Lq = self.arrival_rate, self.calculate_average_queue_length()
        Wq = Lq / arrival_rate
        return Wq

    def calculate_probability_waiting_time_greater_than_t(self, t):
        """
        Calculates the probability that a customer's waiting time exceeds 't'.

        :param t: Time threshold.
        """
        arrival_rate, service_rate, servers, p0 = (
            self.arrival_rate, self.service_rate, self.servers, self.calculate_probability_zero_customers()
        )
        pT = math.exp((-1) * service_rate * t * (
                1
                + ((arrival_rate / service_rate) ** servers) * p0 * (
                    1 - math.exp((-1) * service_rate * t * (servers - arrival_rate / service_rate))) / (
                        math.factorial(servers) * (1 - arrival_rate / (servers * service_rate)) * (
                            servers - 1 - arrival_rate))))
        return pT


class MG1:
    def __init__(self, arrival_rate, service_rate, variance):
        """
        Initializes an M/G/1 queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of the server.
        :param variance: Variance of service times.
        """
        self.arrival_rate = float(arrival_rate)
        self.service_rate = float(service_rate)
        self.variance = float(variance)
        self.name = "MG1"

    def calculate_traffic_intensity(self):
        """Calculates the traffic intensity (ro) of the queueing system."""
        ro = self.arrival_rate / self.service_rate
        return ro

    def calculate_probability_zero_customers(self):
        """Calculates the probability of having zero customers in the system (P0)."""
        ro = self.calculate_traffic_intensity()
        p0 = 1 - ro
        return p0

    def calculate_average_queue_length(self):
        """Calculates the average number of customers in the queue (Lq)."""
        arrival_rate, variance, p0 = self.arrival_rate, self.variance, self.calculate_probability_zero_customers()
        Lq = ((arrival_rate ** 2) * variance + (1 - p0) ** 2) / (2 * p0)
        return Lq

    def calculate_average_system_length(self):
        """Calculates the average number of customers in the system (L)."""
        Lq, ro = self.calculate_average_queue_length(), self.calculate_traffic_intensity()
        L = Lq + ro
        return L

    def calculate_queue_time(self):
        """Calculates the average time a customer spends in the queue (Wq)."""
        Lq, arrival_rate = self.calculate_average_queue_length(), self.arrival_rate
        Wq = Lq / arrival_rate
        return Wq

    def calculate_system_time(self):
        """Calculates the average time a customer spends in the system (W)."""
        Wq, service_rate = self.calculate_queue_time(), self.service_rate
        W = Wq + 1 / service_rate
        return W


class MD1(MG1):
    def __init__(self, arrival_rate, service_rate):
        """
        Initializes an M/D/1 queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of the server.
        """
        super().__init__(arrival_rate, service_rate, 0)  # Variance for M/D/1 is 0
        self.name = "MD1"


class MEk1(MG1):
    def __init__(self, arrival_rate, service_rate, k):
        """
        Initializes an M/Ek/1 queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of the server.
        :param k: Number of phases in Erlang distribution.
        """
        super().__init__(arrival_rate, service_rate, 1 / (k * (service_rate ** 2)))
        self.k = float(k)
        self.name = "MEk1"


class MMSN:
    def __init__(self, arrival_rate, service_rate, servers, N):
        """
        Initializes an M/M/S/N queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of each server.
        :param servers: Number of servers.
        :param N: System capacity.
        """
        self.arrival_rate = float(arrival_rate)
        self.service_rate = float(service_rate)
        self.servers = int(servers)
        self.N = N
        self.name = "MMSN"

    def calculate_traffic_intensity(self):
        """Calculates the traffic intensity (ro) of the queueing system."""
        l, m, s = self.arrival_rate, self.service_rate, self.servers
        ro = l / (m * s)
        return ro

    def calculate_probability_zero_customers(self):
        """Calculates the probability of having zero customers in the system (P0)."""
        somatorio1 = sum((1.0 / math.factorial(i)) * ((self.arrival_rate / self.service_rate) ** i) *
                         (math.factorial(self.N) / math.factorial(self.N - i)) for i in range(self.servers))
        somatorio2 = sum((math.factorial(self.N) * ((self.arrival_rate / self.service_rate) ** i)) /
                         (math.factorial(self.N - i) * math.factorial(self.servers) * (self.servers ** (i - self.servers)))
                         for i in range(self.servers, self.N + 1))
        r = 1.0 / (somatorio1 + somatorio2)
        return r

    def calculate_probability_n_customers(self, n):
        """
        Calculates the probability of having 'n' customers in the system.

        :param n: Number of customers.
        """
        n = int(n)
        if n <= self.servers:
            r = (math.factorial(self.N) * ((self.arrival_rate / self.service_rate) ** n) * self.calculate_probability_zero_customers()) / (
                    math.factorial(self.N - n) * math.factorial(n))
        else:
            if n > self.N:
                r = 0
            else:
                r = (math.factorial(self.N) * ((self.arrival_rate / self.service_rate) ** n) * self.calculate_probability_zero_customers()) / (
                        math.factorial(self.N - n) * math.factorial(self.servers) * (self.servers ** (n - self.servers)))
        return r

    def calculate_queue_probability(self):
        """Calculates the probability of a customer being in the queue (P_queue)."""
        s = self.servers
        S = sum(self.calculate_probability_n_customers(i) for i in range(s))
        return 1 - S

    def calculate_average_queue_length(self):
        """Calculates the average number of customers in the queue (Lq)."""
        r = sum((i - self.servers) * self.calculate_probability_n_customers(i) for i in range(self.servers, self.N + 1))
        return r

    def calculate_average_system_length(self):
        """Calculates the average number of customers in the system (L)."""
        r1 = sum((self.servers - i) * self.calculate_probability_n_customers(i) for i in range(self.servers + 1))
        yeq = self.service_rate * (self.servers - r1)
        r = self.calculate_average_queue_length() + yeq / self.service_rate
        return r

    def calculate_system_time(self):
        """Calculates the average time a customer spends in the system (W)."""
        r1 = sum((self.servers - i) * self.calculate_probability_n_customers(i) for i in range(self.servers + 1))
        yeq = self.service_rate * (self.servers - r1)
        r = self.calculate_average_system_length() / yeq
        return r

    def calculate_queue_time(self):
        """Calculates the average time a customer spends in the queue (Wq)."""
        r1 = sum((self.servers - i) * self.calculate_probability_n_customers(i) for i in range(self.servers + 1))
        yeq = self.service_rate * (self.servers - r1)
        r = self.calculate_average_queue_length() / yeq
        return r


class MMSK:
    def __init__(self, arrival_rate, service_rate, servers, capacity):
        """
        Initializes an M/M/S/K queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of each server.
        :param servers: Number of servers.
        :param capacity: System capacity.
        """
        self.arrival_rate = float(arrival_rate)
        self.service_rate = float(service_rate)
        self.servers = int(servers)
        self.capacity = int(capacity)
        self.name = "MMSK"

    def calculate_traffic_intensity(self):
        """
        Calculates the traffic intensity (ro) of the queueing system.
        """
        l, m, s = self.arrival_rate, self.service_rate, self.servers
        ro = l / (m * s)
        return ro

    def calculate_probability_zero_customers(self):
        """
        Calculates the probability of having zero customers in the system (P0).
        """
        somatorio1 = sum((1.0 / math.factorial(i)) * ((self.arrival_rate / self.service_rate) ** i) *
                         (math.factorial(self.capacity) / math.factorial(self.capacity - i)) for i in range(self.servers))
        somatorio2 = sum((math.factorial(self.capacity) * ((self.arrival_rate / self.service_rate) ** i)) /
                         (math.factorial(self.capacity - i) * math.factorial(self.servers) * (self.servers ** (i - self.servers)))
                         for i in range(self.servers, self.capacity + 1))
        r = 1.0 / (somatorio1 + somatorio2)
        return r

    def calculate_probability(self, n):
        """
        Calculates the probability of having 'n' customers in the system (Pn).

        :param n: Number of customers.
        """
        n = int(n)
        if n <= self.servers:
            r = (1.0 / math.factorial(n)) * ((self.arrival_rate / self.service_rate) ** n) * self.calculate_probability_zero_customers()
        else:
            if n > self.capacity:
                r = 0
            else:
                r = (1.0 / (math.factorial(self.servers) * (self.servers ** (n - self.servers)))) * (
                        (self.arrival_rate / self.service_rate) ** n) * self.calculate_probability_zero_customers()
        return r

    def calculate_queue_probability(self):
        """
        Calculates the probability of the system being in the queue state (Pqueue).
        """
        s = self.servers
        S = sum(self.calculate_probability(i) for i in range(s))
        queue_probability = 1 - S
        return queue_probability

    def calculate_expected_customers_queue(self):
        """
        Calculates the expected number of customers in the queue (Nqueue).
        """
        r1 = ((self.arrival_rate / self.service_rate) ** self.servers)
        r2 = (self.arrival_rate / (self.service_rate * self.servers))
        r3 = math.factorial(self.servers) * ((1 - r2) ** (2))
        r4 = 1.0 - (r2 ** (self.capacity - self.servers)) - (self.capacity - self.servers) * (
                r2 ** (self.capacity - self.servers)) * (1 - r2)
        expected_customers_queue = (r1 * r2 * r4 * self.calculate_probability_zero_customers()) / r3
        return expected_customers_queue

    def calculate_expected_customers_system(self):
        """
        Calculates the expected number of customers in the entire system (Nsystem).
        """
        r1 = sum(i * self.calculate_probability(i) for i in range(0, self.servers))
        r2 = sum(self.calculate_probability(i) for i in range(0, self.servers))
        expected_customers_system = r1 + self.calculate_expected_customers_queue() + (self.servers * (1 - r2))
        return expected_customers_system

    def calculate_system_time(self):
        """
        Calculates the average time a customer spends in the system (Tsystem).
        """
        system_time = self.calculate_expected_customers_system() / (self.arrival_rate * (1 - self.calculate_probability(self.capacity)))
        return system_time

    def calculate_queue_time(self):
        """
        Calculates the average time a customer spends in the queue (Tqueue).
        """
        queue_time = self.calculate_expected_customers_queue() / (self.arrival_rate * (1 - self.calculate_probability(self.capacity)))
        return queue_time


class MG1K:
    def __init__(self, arrival_rate, service_rate, variance, capacity):
        """
        Initializes an M/G/1/K queueing system.

        :param arrival_rate: Arrival rate to the system.
        :param service_rate: Service rate of the server.
        :param variance: Variance of service times.
        :param capacity: System capacity.
        """
        self.arrival_rate = float(arrival_rate)
        self.service_rate = float(service_rate)
        self.expectation = 1 / service_rate
        self.variance = variance
        self.capacity = capacity
        self.name = "MG1"

    def calculate_traffic_intensity(self):
        """
        Calculates the traffic intensity (ro) of the queueing system.
        """
        ro = self.arrival_rate / self.service_rate
        return ro

    def calculate_probability_k(self):
        """
        Calculates the probability of having 'K' or more customers in the system (Pk).
        """
        s2, k, ro = self.variance, self.capacity, self.calculate_traffic_intensity()
        pk = (ro ** ((2 + math.sqrt(ro) * (s2) - math.sqrt(ro) + 2 * (k)) /
                     (2 + math.sqrt(ro) * (s2) - math.sqrt(ro)))) * (ro - 1)
        pk = pk / (ro ** (
                2 * (2 + math.sqrt(ro) * (s2) - math.sqrt(ro) + (k)) / (2 + math.sqrt(ro) * (s2) - math.sqrt(ro)) - 1))
        return pk

    def calculate_throughput(self):
        """
        Calculates the system throughput (O).
        """
        l, pk = self.arrival_rate, self.calculate_probability_k()
        throughput = l * (1 - pk)
        return throughput

    def calculate_mi_h(self):
        """
        Calculates the effective service rate for high traffic intensity (Mi_h).
        """
        mi, s2 = self.service_rate, self.variance
        mi_h = 2 * mi / (1 + s2)  # to testando outro mi_h 2*mi/(1+s2*(mi**2))
        return mi_h

    def calculate_mi_l_h(self):
        """
        Calculates the effective service rate for low to high traffic intensity transition (Mi_l_h).
        """
        pk_l, mi_h = self.calculate_probability_k(), self.calculate_mi_h()
        mi_l_h = (1 - pk_l) * mi_h
        return mi_l_h

    def calculate_probability_k_l(self):
        """
        Calculates the probability of having 'K' or more customers in the system for low to high traffic intensity transition (Pk_l).
        """
        lamb, mi, k, mi_h, r1, r2 = self.arrival_rate, self.service_rate, self.capacity, self.calculate_mi_h(), *self.calculate_r()
        pk_h = (((mi + mi_h) / mi_h) -
                (lamb * ((r2 ** k - r1 ** k) - (r2 ** (k - 1) - r1 ** (k - 1)))) /
                (mi_h * ((r2 ** (k + 1) - r1 ** (k + 1)) - (r2 ** k - r1 ** k)))) ** (-1)
        return pk_h

    def calculate_r(self):
        """
        Calculates the values of r1 and r2 for low to high traffic intensity transition.
        """
        lamb, mi_h, z = self.arrival_rate, self.calculate_mi_h(), self.calculate_z()
        r = (lamb + 2 * mi_h) / (2 * mi_h)
        zeta = math.sqrt(z) / (2 * mi_h)
        return r - zeta, r + zeta

    def calculate_z(self):
        """
        Calculates the value of z for low to high traffic intensity transition.
        """
        lamb, mi_h = self.arrival_rate, self.calculate_mi_h()
        z = ((lamb + 2 * mi_h) ** 2) - 4 * lamb * mi_h
        return z
