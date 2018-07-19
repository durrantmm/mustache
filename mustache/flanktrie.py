# Adapted from code by Shubhadeep Roychowdhury available at https://towardsdatascience.com/implementing-a-trie-data-structure-in-python-in-less-than-100-lines-of-code-a877ea23c1a1
from typing import Tuple
import sys

class TrieNode(object):
    """
    Our trie node implementation. Very basic. but does the job
    """

    def __init__(self, char: str, parent=None):
        self.char = char
        self.children = {}
        self.parent = parent
        # Is it the last character of the word.`
        self.word_count = 0
        # How many times this character appeared in the addition process
        self.counter = 1
        self.total_lifetime_children = 0
        self.qual = 0


class Trie:

    def __init__(self):
        self.root = TrieNode('')
        self.total_words = 0
        self.total_nodes = 0
        self.total_base_quality = 0

    def add(self, word: str, quals: list):
        """
        Adding a word in the trie structure
        """

        if len(word) != len(quals):
            print("Added word and character quality vectors ahve different lengths")
            sys.exit()

        node = self.root
        for i in range(len(word)):
            char = word[i]
            qual = quals[i]
            found_in_child = False
            # Search for the character in the children of the present `node`

            if char in node.children:
                child = node.children[char]
                child.counter += 1
                child.qual += qual
                node = child
            else:
                new_node = TrieNode(char, node)
                new_node.qual += qual
                node.children[char] = new_node
                node.total_lifetime_children += 1
                node = new_node

        node.word_count += 1


    def make_subtrie(self, words):
        subtrie = Trie()

        for word in words:
            trie_node = self.root
            subtrie_node = subtrie.root

            for char in word:

                trie_child = trie_node.children[char]

                if char in subtrie_node.children:
                    subtrie_child = subtrie_node.children[char]

                else:
                    subtrie_child = TrieNode(char, subtrie_node)

                    subtrie_child.qual = trie_child.qual
                    subtrie_child.word_count = trie_child.word_count
                    subtrie.total_words += trie_child.word_count

                    subtrie_node.children[char] = subtrie_child
                    subtrie_node.total_lifetime_children = len(subtrie_node.children)


                trie_node = trie_child
                subtrie_node = subtrie_child


        return subtrie


    def find_prefix(self, prefix: str) -> Tuple[bool, int]:
        """
        Check and return
          1. If the prefix exists in any of the words we added so far
          2. If yes then how may words actually have the prefix
        """
        node = self.root
        # If the root node has no children, then return False.
        # Because it means we are trying to search in an empty trie
        if not self.root.children:
            return False, 0
        for char in prefix:
            char_not_found = True
            # Search through all the children of the present `node`
            for child in node.children:
                if child.char == char:
                    # We found the char existing in the child.
                    char_not_found = False
                    # Assign node as the child containing the char and break
                    node = child
                    break
            # Return False anyway when we did not find a char.
            if char_not_found:
                return False, 0
        # Well, we are here means we have found the prefix. Return true to indicate that
        # And also the counter of the last node. This indicates how many words have this
        # prefix
        return True, node.counter

    def traverse_seqs(self, prefix='', node=None):
        if node is None:
            node = self.root

        word = prefix + node.char

        if len(node.children) == 0:
            return [word]

        words = []
        for child in node.children.values():
            words = words + self.traverse_seqs(word, child)

        return words

    def traverse_quals(self, prefix=[], node=None):
        if node is None:
            node = self.root
            qual = []
        else:
            qual = prefix + [node.qual]

        if len(node.children) == 0:
            return [qual]

        quals = []
        for child in node.children.values():
            quals = quals + self.traverse_quals(qual, child)

        return quals

    def traverse_both(self):

        seqs = self.traverse_seqs()
        quals = self.traverse_quals()

        both = list(zip(seqs, quals))
        return(both)

    def calc_total_words(self, word):
        total_words = 0
        node = self.root

        for c in word:
            total_words += node.children[c].word_count
            node = node.children[c]

        return total_words

    def calc_total_words_before_lifetime_child(self, word):
        total_words = 0
        node = self.root

        for c in word:
            if node.total_lifetime_children > 1:
                break
            total_words += node.children[c].word_count
            node = node.children[c]

        return total_words

    def calc_total_shared_words(self, word1, word2):

        shared_words = 0
        node = self.root

        for i in range(min([len(word1), len(word2)])):
            c1, c2 = word1[i], word2[i]
            if c1 != c2:
                break
            else:
                shared_words += node.children[c1].word_count

            node = node.children[c1]

        return shared_words

    def calc_total_unique_shared_words(self, word1, word2):

        shared_words = 0
        node = self.root

        for i in range(min([len(word1), len(word2)])):
            c1, c2 = word1[i], word2[i]
            if c1 != c2:
                break
            elif len(node.children) > 1:
                shared_words = node.children[c1].word_count
            else:
                shared_words += node.children[c1].word_count

            node = node.children[c1]

        return shared_words

    def calc_total_unique_words(self, word1, word2):

        unique_words = 0
        node = self.root

        same_path = True
        for i in range(len(word1)):
            c1 = word1[i]
            if i >= len(word2) or c1 != word2[i]:
                same_path = False

            if not same_path:
                unique_words += node.children[c1].word_count

            node = node.children[c1]
        return unique_words

    def delete_word(self, word):
        node = self.root
        for char in word:
            node = node.children[char]

        while True:
            parent_node = node.parent

            if parent_node.char == '' or len(parent_node.children) > 1:
                del parent_node.children[char]
                break

            del parent_node.children[char]
            node = parent_node
            char = node.char

    def calc_word_count_diff(self, word1, word2):

        node = self.root
        #for c in word1:
        #    for


if __name__ == "__main__":
    trie = Trie()
    trie.add('hellothere', [40]*len('hellothere'))
    trie.add('hellothen', [40]*len('hellothen'))
    trie.add('hellohowareyou', [40]*len('hellohowareyou'))
    trie.add('hello', [40]*len('hello'))

    print(trie.calc_total_shared_words('hellothere', 'hellothen'))
    print(trie.calc_total_unique_words('hellothere', 'hellothen'))
    print(trie.calc_total_unique_words('hellothen', 'hellothere'))
    print(trie.calc_total_words_before_lifetime_child('hellothen'))