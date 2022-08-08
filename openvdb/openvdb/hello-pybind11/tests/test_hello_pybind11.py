# This kind of import is automatically done when importing hello from outside
import hello
import unittest


class TestHello(unittest.TestCase):
    def test_hello(self):
        hello.hello()

    def test_return_two(self):
        self.assertEqual(hello.return_two(), 2)

    def test_pass_coord(self):
        coord = (1, 2, 3)
        result = hello.pass_coord(coord)
        self.assertEqual(coord, result)

    def test_pass_coord_default(self):
        default_coord = (0, 0, 0)
        result = hello.pass_coord()
        self.assertEqual(default_coord, result)


if __name__ == "__main__":
    unittest.main()
    # You can run all python test with:
    # ctest -R python -V
    # from the build folder
