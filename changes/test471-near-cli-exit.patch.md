**Headline:** The test for `planetgen query near` with no distance now checks the exit code and the usage message the way the command-line parser reports them.

### Changed

- TEST.471: `test_the_command_line_asks_for_the_nearest_few` expected the message in the exit value; the parser exits with code 2 and writes the reason ("give --distance ...") to stderr, so the test failed on main. It now asserts the exit code and reads the message from captured stderr, as TEST.465 does. Test only.
