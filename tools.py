"""
Convenience functions and helpers.
"""

import cPickle

import numpy as np


class Result:
    """
    Result class.
    """
    RESULT_TYPES = ['ph_ph_mi', 'ph_amp_mi', 'ph_ph_caus', 'ph_amp_caus']

    def __init__(self, array, metric, algorithm):
        assert metric in self.RESULT_TYPES, 'Wrong metric: %s' % metric
        assert isinstance(array, np.ndarray)

        self.result = array
        self.metric = metric
        self.algorithm = algorithm

    @property
    def key(self):
        return "%s_%s" % (self.metric, self.algorithm)

    def __str__(self):
        return "Contains %s result with %s algorithm of shape: %s" % (
            self.metric, self.algorithm, str(self.result.shape))


class ResultsContainer(Result):
    """
    Collection of Results.
    """
    @staticmethod
    def validate_single_result(result):
        """
        Validate single result.
        """
        assert isinstance(result, (list, tuple)), (
            'Wrong type: %s' % type(result))
        assert len(result) == 4, ('Wrong length: %s' % len(result))
        for single_result in result:
            assert isinstance(single_result, dict)
            assert len(list(single_result.keys())) == 2
            for key, value in single_result.iteritems():
                assert isinstance(value, np.ndarray)

    @classmethod
    def from_algorithm_output(cls, results_list):
        cls.validate_single_result(results_list)
        output = []
        for result, metric in zip(results_list, cls.RESULT_TYPES):
            for algorithm in result:
                output.append(Result(
                    result[algorithm],
                    metric=metric,
                    algorithm=algorithm))

        return cls(output)

    @classmethod
    def from_saved_file(cls, filename):
        print('Loading from %s...' % filename)
        with open(filename, "rb") as f:
            raw_data = cPickle.load(f)

        def get_metric_algorithm(key):
            for metric in cls.RESULT_TYPES:
                if metric in key:
                    algorithm = key.replace("%s_" % metric, "")
                    return metric, algorithm

        results = []
        for key, array in raw_data.iteritems():
            metric, algorithm = get_metric_algorithm(key)
            results.append(Result(array, metric=metric, algorithm=algorithm))
        return cls(results)

    def __init__(self, results):
        assert all(isinstance(result, Result) for result in results)
        self.results = results

    def __getitem__(self, key):
        metric, algo = key
        for result in self.results:
            if result.metric == metric and result.algorithm == algo:
                return result.result

    def save(self, filename):
        saving_dict = {
            result.key: result.result for result in self.results
        }
        with open(filename, "wb") as f:
            cPickle.dump(saving_dict, f, protocol=cPickle.HIGHEST_PROTOCOL)
        print('Saving done to %s' % filename)


class SurrogatesContainer(ResultsContainer):
    """
    Results from surrogates.
    """
    @classmethod
    def from_algorithm_output(cls, results_list):
        output_dict = {}
        for single_result in results_list:
            cls.validate_single_result(single_result)
            for result, metric in zip(single_result, cls.RESULT_TYPES):
                for algorithm in result:
                    key = "%s_%s" % (metric, algorithm)
                    # if we have the key, stack results
                    if key in output_dict:
                        output_dict[key] = np.dstack([output_dict[key],
                                                      result[algorithm]])
                    # if not create new key
                    else:
                        output_dict[key] = result[algorithm]

        output = []
        for key, array in output_dict.iteritems():
            metric, algorithm = key.split("_")
            output.append(Result(array, metric=metric, algorithm=algorithm))

    return cls(output)
