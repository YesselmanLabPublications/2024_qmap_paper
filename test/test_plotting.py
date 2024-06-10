import pytest

from qmap_paper.plotting import colors_for_sequence, find_stretches, trim


class TestColorsForSequence:
    def test_valid_sequence(self):
        """
        Test valid DNA sequences with expected color mappings.
        """
        assert colors_for_sequence("A") == ["red"]
        assert colors_for_sequence("C") == ["blue"]
        assert colors_for_sequence("G") == ["orange"]
        assert colors_for_sequence("T") == ["green"]
        assert colors_for_sequence("ACGT") == ["red", "blue", "orange", "green"]
        assert colors_for_sequence("TGCA") == ["green", "orange", "blue", "red"]

    def test_empty_sequence(self):
        """
        Test an empty DNA sequence.
        """
        assert colors_for_sequence("") == []

    def test_invalid_character_in_sequence(self):
        """
        Test that invalid characters in the DNA sequence raise a ValueError.
        """
        with pytest.raises(ValueError):
            colors_for_sequence("X")
        with pytest.raises(ValueError):
            colors_for_sequence("A1C")
        with pytest.raises(ValueError):
            colors_for_sequence("AGCTZ")

    def test_mixed_valid_and_invalid_sequence(self):
        """
        Test sequences with a mix of valid and invalid characters.
        """
        with pytest.raises(ValueError):
            colors_for_sequence("AGXT")
        with pytest.raises(ValueError):
            colors_for_sequence("XYZ")

    def test_repetitive_valid_sequence(self):
        """
        Test sequences with repetitive valid characters.
        """
        assert colors_for_sequence("AAAA") == ["red", "red", "red", "red"]
        assert colors_for_sequence("CCCC") == ["blue", "blue", "blue", "blue"]
        assert colors_for_sequence("GGGG") == ["orange", "orange", "orange", "orange"]
        assert colors_for_sequence("TTTT") == ["green", "green", "green", "green"]


class TestFindStretches:
    """Test cases for the find_stretches function."""

    def test_regular_cases(self):
        """Test cases with regular inputs to check correct identification of consecutive stretches."""
        assert find_stretches([3, 4, 5, 10, 11, 12]) == [[3, 5], [10, 12]]
        assert find_stretches([1, 2, 3, 7, 8, 10]) == [[1, 3], [7, 8], [10, 10]]

    def test_single_stretch(self):
        """Test cases where the list contains only one stretch of consecutive numbers."""
        assert find_stretches([4, 5, 6, 7]) == [[4, 7]]

    def test_no_consecutive(self):
        """Test cases where there are no consecutive numbers in the list."""
        assert find_stretches([1, 3, 5, 7]) == [[1, 1], [3, 3], [5, 5], [7, 7]]

    def test_empty_list(self):
        """Test the function with an empty list."""
        assert find_stretches([]) == []

    def test_single_element_list(self):
        """Test case with a single element in the list."""
        assert find_stretches([42]) == [[42, 42]]

    def test_descending_order(self):
        """Test case where the elements are in descending order."""
        assert find_stretches([20, 19, 18, 17, 16]) == [[16, 20]]

    def test_duplicates(self):
        """Test case where the list contains duplicate elements."""
        assert find_stretches([3, 3, 4, 5, 6, 6, 7]) == [[3, 7]]

    def test_unordered_elements(self):
        """Test case with elements in random order."""
        assert find_stretches([10, 1, 3, 2, 5, 4, 6, 8, 9, 7]) == [[1, 10]]


class TestTrimFunction:

    def test_trim_string(self):
        assert trim("hello world", 1, 1) == "ello worl"
        assert trim("abcdefg", 2, 2) == "cde"
        assert trim("12345", 0, 0) == "12345"

    def test_trim_list(self):
        assert trim([1, 2, 3, 4, 5], 1, 1) == [2, 3, 4]
        assert trim(["a", "b", "c", "d"], 2, 1) == ["c"]
        assert trim([10, 20, 30, 40, 50], 0, 0) == [10, 20, 30, 40, 50]

    def test_invalid_content_type(self):
        with pytest.raises(
            TypeError, match="Invalid content type. Please provide a string or a list."
        ):
            trim(12345, 1, 1)
        with pytest.raises(
            TypeError, match="Invalid content type. Please provide a string or a list."
        ):
            trim({1: "a", 2: "b"}, 1, 1)

    def test_boundary_cases(self):
        assert trim("", 1, 1) == ""
        assert trim([], 1, 1) == []
        assert trim("a", 1, 1) == ""
        assert trim(["a"], 1, 1) == []
