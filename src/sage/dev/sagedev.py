- David Roe, Frej Drejhammar, Julian Rueth, Martin Raum, Nicolas M. Thiery, R.,
#                          Volker Braun <vbraun.name@gmail.com>
import os
import urllib, urlparse
import re

from patch import MercurialPatchMixin
from sage.env import TRAC_SERVER_URI
# http://stackoverflow.com/questions/12093748/how-do-i-check-for-valid-git-branch-names
GIT_BRANCH_REGEX = re.compile(
    r'^(?!.*/\.)(?!.*\.\.)(?!/)(?!.*//)(?!.*@\{)(?!.*\\)'
    r'[^\040\177 ~^:?*[]+(?<!\.lock)(?<!/)(?<!\.)$')
#
# The first line should contain a short summary of your changes, the
# following lines should contain a more detailed description. Lines
# starting with '#' are ignored.
class SageDev(MercurialPatchMixin):
    - ``config`` -- a :class:`~sage.dev.config.Config` or ``None``
      (default: ``None``), the configuration of this object; the
      defaults uses the configuration stored in ``DOT_SAGE/devrc``.
    - ``UI`` -- a :class:`~sage.dev.user_interface.UserInterface` or ``None`` (default:
            self.git = GitInterface(self.config['git'], self._UI)
                self._UI.show('The developer scripts used to store some of their data in "{0}".'
                              ' This file has now moved to "{1}". I moved "{0}" to "{1}". This might'
                              ' cause trouble if this is a fresh clone of the repository in which'
                              ' you never used the developer scripts before. In that case you'
                              ' should manually delete "{1}" now.', old_file, new_file)
        move_legacy_saving_dict('ticketfile', self.config['sagedev'].get(
            'ticketfile', os.path.join(DOT_SAGE, 'branch_to_ticket')), ticket_file)
        move_legacy_saving_dict('branchfile', self.config['sagedev'].get(
            'branchfile', os.path.join(DOT_SAGE, 'ticket_to_branch')), branch_file)
        move_legacy_saving_dict('dependenciesfile', self.config['sagedev'].get(
            'dependenciesfile', os.path.join(DOT_SAGE, 'dependencies')), dependencies_file)
        move_legacy_saving_dict('remotebranchesfile', self.config['sagedev'].get(
            'remotebranchesfile', os.path.join(DOT_SAGE, 'remote_branches')), remote_branches_file)
            from sage.dev.misc import tmp_dir
            self._tmp_dir = tmp_dir()
    def create_ticket(self):
        Create a new ticket on trac.
            :meth:`checkout`, :meth:`pull`, :meth:`edit_ticket`
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)

            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
        This fails if the internet connection is broken::
            sage: dev.trac._connected = False
            sage: UI.append("Summary: ticket7\ndescription")
            sage: dev.create_ticket()
            A network error ocurred, ticket creation aborted.
            Your command failed because no connection to trac could be established.
            sage: dev.trac._connected = True
        """
        try:
            ticket = self.trac.create_ticket_interactive()
        except OperationCancelledError:
            self._UI.debug("Ticket creation aborted.")
            raise
        except TracConnectionError as e:
            self._UI.error("A network error ocurred, ticket creation aborted.")
            raise
        ticket_url = urlparse.urljoin(self.trac._config.get('server', TRAC_SERVER_URI), str(ticket))
        self._UI.show("Created ticket #{0} at {1}.".format(ticket, ticket_url))
        self._UI.info(['',
                       '(use "{0}" to create a new local branch)'
                       .format(self._format_command("checkout", ticket=ticket))])
        return ticket
    def checkout(self, ticket=None, branch=None, base=''):
        r"""
        Checkout another branch.
        If ``ticket`` is specified, and ``branch`` is an existing local branch,
        then ``ticket`` will be associated to it, and ``branch`` will be
        checked out into the working directory.
        Otherwise, if there is no local branch for ``ticket``, the branch
        specified on trac will be pulled to ``branch`` unless ``base`` is
        set to something other than the empty string ``''``. If the trac ticket
        does not specify a branch yet or if ``base`` is not the empty string,
        then a new one will be created from ``base`` (per default, the master
        branch).
        If ``ticket`` is not specified, then checkout the local branch
        ``branch`` into the working directory.
        INPUT:
        - ``ticket`` -- a string or an integer identifying a ticket or ``None``
          (default: ``None``)
        - ``branch`` -- a string, the name of a local branch; if ``ticket`` is
          specified, then this defaults to ticket/``ticket``.
        - ``base`` -- a string or ``None``, a branch on which to base a new
          branch if one is going to be created (default: the empty string
          ``''`` to create the new branch from the master branch), or a ticket;
          if ``base`` is set to ``None``, then the current ticket is used. If
          ``base`` is a ticket, then the corresponding dependency will be
          added. Must be ``''`` if ``ticket`` is not specified.
        .. SEEALSO::
            :meth:`pull`, :meth:`create_ticket`, :meth:`vanilla`
        TESTS:
        Set up a single user for doctesting::
            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
        Create a few branches::
            sage: dev.git.silent.branch("branch1")
            sage: dev.git.silent.branch("branch2")
        Checking out a branch::
            sage: dev.checkout(branch="branch1")
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.current_branch()
            'branch1'
        Create a ticket and checkout a branch for it::
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.current_branch()
            'ticket/1'
        """
        if ticket is not None:
            self.checkout_ticket(ticket=ticket, branch=branch, base=base)
        elif branch is not None:
            if base != '':
                raise SageDevValueError("base must not be specified if no ticket is specified.")
            self.checkout_branch(branch=branch)
        else:
            raise SageDevValueError("at least one of ticket or branch must be specified.")

        ticket = self._current_ticket()
        branch = self.git.current_branch()
        if ticket:
            self._UI.show(['On ticket #{0} with associated local branch "{1}".'], ticket, branch)
        else:
            self._UI.show(['On local branch "{0}" without associated ticket.'], branch)
        self._UI.info(['', 
                       'Use "{0}" to save changes in a new commit when you are finished editing.'],
                      self._format_command("commit"))
    def checkout_ticket(self, ticket, branch=None, base=''):
        Checkout the branch associated to ``ticket``.
        associated to it, and ``branch`` will be checked out into the working directory.
        specified on trac will be pulled to ``branch`` unless ``base`` is
            :meth:`pull`, :meth:`create_ticket`, :meth:`vanilla`
        Alice tries to checkout ticket #1 which does not exist yet::
            sage: alice.checkout(ticket=1)
            Ticket name "1" is not valid or ticket does not exist on trac.
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Now alice can check it out, even though there is no branch on the
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        If Bob commits something to the ticket, a ``checkout`` by Alice
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        If Alice had not checked that ticket out before, she would of course
            sage: alice.checkout(ticket=1) # ticket #1 refers to the non-existant branch 'ticket/1'
            Ticket #1 refers to the non-existant local branch "ticket/1". If you have not
            manually interacted with git, then this is a bug in sagedev. Removing the
            association from ticket #1 to branch "ticket/1".
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Checking out a ticket with untracked files::
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Checking out a ticket with untracked files which make a checkout
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.checkout(ticket=1)
            This happened while executing "git -c user.email=doc@test.test -c
            user.name=alice checkout ticket/1".
        Checking out a ticket with uncommited changes::
            sage: open("tracked", "w").close()
            sage: alice.checkout(ticket=2)
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Keep/stash] d
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.

        Now follow some single user tests to check that the parameters are interpreted correctly::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("_dependencies_for_ticket")

        First, create some tickets::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: UI.append("Summary: ticket2\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.silent.commit(allow_empty=True, message="second commit")
            sage: dev.git.commit_for_branch('ticket/2') != dev.git.commit_for_branch('ticket/1')
            True

        Check that ``base`` works::

            sage: UI.append("Summary: ticket3\ndescription")
            sage: dev.create_ticket()
            Created ticket #3 at https://trac.sagemath.org/3.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=3" to create a new local branch)
            3
            sage: dev.checkout(ticket=3, base=2)
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.commit_for_branch('ticket/3') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(3)
            (2,)
            sage: UI.append("Summary: ticket4\ndescription")
            sage: dev.create_ticket()
            Created ticket #4 at https://trac.sagemath.org/4.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=4" to create a new local branch)
            4
            sage: dev.checkout(ticket=4, base='ticket/2')
            On ticket #4 with associated local branch "ticket/4".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.commit_for_branch('ticket/4') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(4)
            ()

        In this example ``base`` does not exist::

            sage: UI.append("Summary: ticket5\ndescription")
            sage: dev.create_ticket()
            Created ticket #5 at https://trac.sagemath.org/5.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=5" to create a new local branch)
            5
            sage: dev.checkout(ticket=5, base=1000)
            Ticket name "1000" is not valid or ticket does not exist on trac.
        In this example ``base`` does not exist locally::

            sage: UI.append("Summary: ticket6\ndescription")
            sage: dev.create_ticket()
            Created ticket #6 at https://trac.sagemath.org/6.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=6" to create a new local branch)
            6
            sage: dev.checkout(ticket=6, base=5)
            Branch field is not set for ticket #5 on trac.

        Creating a ticket when in detached HEAD state::

            sage: dev.git.super_silent.checkout('HEAD', detach=True)
            sage: UI.append("Summary: ticket detached\ndescription")
            sage: dev.create_ticket()
            Created ticket #7 at https://trac.sagemath.org/7.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=7" to create a new local branch)
            7
            sage: dev.checkout(ticket=7)
            On ticket #7 with associated local branch "ticket/7".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.current_branch()
            'ticket/7'

        Creating a ticket when in the middle of a merge::

            sage: dev.git.super_silent.checkout('-b','merge_branch')
            sage: with open('merge', 'w') as f: f.write("version 0")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','some change')
            sage: dev.git.super_silent.checkout('ticket/7')
            sage: with open('merge', 'w') as f: f.write("version 1")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','conflicting change')
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.silent.merge('merge_branch')
            ....: except GitError: pass
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #8 at https://trac.sagemath.org/8.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=8" to create a new local branch)
            8
            sage: UI.append("cancel")
            sage: dev.checkout(ticket=8)
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] cancel
            Aborting checkout of branch "ticket/8".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)
            sage: dev.git.reset_to_clean_state()

        Creating a ticket with uncommitted changes::

            sage: open('tracked', 'w').close()
            sage: dev.git.silent.add('tracked')
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #9 at https://trac.sagemath.org/9.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=9" to create a new local branch)
            9

        The new branch is based on master which is not the same commit
        as the current branch ``ticket/7``, so it is not a valid
        option to ``'keep'`` changes::

            sage: UI.append("cancel")
            sage: dev.checkout(ticket=9)
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Aborting checkout of branch "ticket/9".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)

        Finally, in this case we can keep changes because the base is
        the same commit as the current branch::

            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #10 at https://trac.sagemath.org/10.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=10" to create a new local branch)
            10
            sage: UI.append("keep")
            sage: dev.checkout(ticket=10, base='ticket/7')
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Keep/stash] keep
            On ticket #10 with associated local branch "ticket/10".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        # if branch points to an existing branch make it the ticket's branch and check it out
            self._UI.debug('The branch for ticket #{0} is now "{1}".', ticket, branch)
            self._UI.debug('Now checking out branch "{0}".', branch)
            self.checkout_branch(branch)
        # if there is a branch for ticket locally, check it out
                self._UI.debug('Checking out branch "{0}".', branch)
                self.checkout_branch(branch)
            raise SageDevValueError('currently on no ticket, "base" must not be None')
            base = self._local_branch_for_ticket(base, pull_if_not_found=True)
                    self._UI.debug('The branch field on ticket #{0} is not set. Creating a new branch'
                                   ' "{1}" off the master branch "{2}".', ticket, branch, MASTER_BRANCH)
                    # pull the branch mentioned on trac
                        self._UI.error('The branch field on ticket #{0} is set to the non-existent "{1}".'
                                       ' Please set the field on trac to a field value.', 
                                       ticket, remote_branch)
                        self._UI.info(['', '(use "{0}" to edit the ticket description)'],
                                       self._format_command("edit-ticket", ticket=ticket))
                        self.pull(remote_branch, branch)
                        self._UI.debug('Created a new branch "{0}" based on "{1}".',
                                       branch, remote_branch)
                        self._UI.error('Could not check out ticket #{0} because the remote branch "{1}"'
                                       ' for that ticket could not be pulled (network connection?).', 
                                       ticket, remote_branch)
                    self._UI.show('About to create a new branch for #{0} based on "{1}". However, the trac'
                                  ' ticket for #{0} already refers to the branch "{2}". The new branch will'
                                  ' not contain any work that has already been done on "{2}".',
                                  ticket, base, remote_branch)
                    if not self._UI.confirm('Create fresh branch?', default=False):
                        command += self._format_command("checkout", ticket=ticket)
                        self._UI.info(['', 'Use "{1}" to work on a local copy of the existing remote branch "{0}".'], 
                                      remote_branch, command)
                self._UI.debug('Creating a new branch for #{0} based on "{1}".', ticket, base)
                self._UI.debug('Deleting local branch "{0}".', branch)
            self._UI.debug("Locally recording dependency on {0} for #{1}.",
                           ", ".join(["#"+str(dep) for dep in dependencies]), ticket)
        self._set_remote_branch_for_branch(branch, self._remote_branch_for_ticket(ticket))
        self._UI.debug('Checking out to newly created branch "{0}".'.format(branch))
        self.checkout_branch(branch)
    def checkout_branch(self, branch, helpful=True):
        Checkout to the local branch ``branch``.
        - ``branch`` -- a string, the name of a local branch
        Checking out a branch::
            sage: dev.checkout(branch="branch1")
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.checkout(branch="branch3")
            Branch "branch3" does not exist locally.
            <BLANKLINE>
            #  (use "sage --dev tickets" to list local branches)
        Checking out branches with untracked files::
            sage: open("untracked", "w").close()
            sage: dev.checkout(branch="branch2")
            On local branch "branch2" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Checking out a branch with uncommitted changes::
            sage: open("tracked", "w").close()
            sage: UI.append("cancel")
            sage: dev.checkout(branch="branch1")
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Aborting checkout of branch "branch1".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)
            sage: dev.checkout(branch="branch1")
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] s
            Your changes have been moved to the git stash stack. To re-apply your changes
            later use "git stash apply".
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.

        And retrieve the stashed changes later::

            sage: dev.checkout(branch='branch2')
            On local branch "branch2" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.git.echo.stash('apply')
            # On branch branch2
            # Changes not staged for commit:
            #   (use "git add <file>..." to update what will be committed)
            #   (use "git checkout -- <file>..." to discard changes in working directory)
            #
            #   modified:   tracked
            #
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            no changes added to commit (use "git add" and/or "git commit -a")
            sage: UI.append("discard")
            sage: dev.checkout(branch="branch1")
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Checking out a branch when in the middle of a merge::
            sage: UI.append('r')
            sage: dev.checkout(branch='merge_branch')
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] r
            On local branch "merge_branch" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.

        Checking out a branch when in a detached HEAD::
            sage: dev.checkout(branch='branch1')
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.checkout(branch='branch1')
            <BLANKLINE>
                tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Checking out a branch with untracked files that would be overwritten by
        the checkout::
            sage: dev.checkout(branch='branch2')
            This happened while executing "git -c user.email=doc@test.test -c
            user.name=doctest checkout branch2".
            error: The following untracked working tree files would be overwritten
            by checkout:
            self.reset_to_clean_state(helpful=False)
            if helpful:
                self._UI.show('Aborting checkout of branch "{0}".', branch)
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
            self.clean(error_unless_clean=(current_commit != target_commit))
            if helpful:
                self._UI.show('Aborting checkout of branch "{0}".', branch)
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
            # this leaves locally modified files intact (we only allow this to happen
            # if current_commit == target_commit
    def pull(self, ticket_or_remote_branch=None, branch=None):
        Pull ``ticket_or_remote_branch`` to ``branch``.
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
        Bob attempts to pull for the ticket but fails because there is no
            sage: bob.pull(1)
            Branch field is not set for ticket #1 on trac.
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
        Bob pulls the changes for ticket 1::
            sage: bob.pull()
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
        Alice can now pull the changes by Bob without the need to merge
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: bob: added alices_file
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
        Now, the pull fails; one would have to use :meth:`merge`::
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort
        Undo the latest commit by alice, so we can pull again::
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: bob: added bobs_other_file
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
                    bobs_other_file
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort
            raise SageDevValueError('No "ticket_or_remote_branch" specified to pull.')
        self._UI.debug('Fetching remote branch "{0}" into "{1}".'.format(remote_branch, branch))
            self.merge(remote_branch, pull=True)
                # then just nothing happened and we can abort the pull
                e.explain = 'Fetching "{0}" into "{1}" failed.'.format(remote_branch, branch)
                    e.explain += " Most probably this happened because the fetch did not" \
                                 " resolve as a fast-forward, i.e., there were conflicting changes."
                    e.advice = 'You can try to use "{2}" to checkout "{1}" and then use "{3}"' \
                               ' to resolve these conflicts manually.'.format(
                                   remote_branch, branch, 
                                   self._format_command("checkout", branch=branch), 
                                   self._format_command("merge", remote_branch, pull=True))
                    e.explain += "We did not expect this case to occur.  If you can explain" \
                                 " your context in sage.dev.sagedev it might be useful to others."
            - :meth:`push` -- Push changes to the remote server.  This
              is the next step once you've committed some changes.
            - :meth:`diff` -- Show changes that will be committed.
            sage: dev.git.super_silent.checkout('-b', 'branch1')
            sage: dev._UI.extend(["y", "added tracked", "y", "y"])
            <BLANKLINE>
                tracked
            <BLANKLINE>
            Start tracking any of these files? [yes/No] y
            Start tracking "tracked"? [yes/No] y
            Commit your changes to branch "branch1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are
            #  done.
            sage: with open("tracked", "w") as F: F.write("foo")
            sage: dev._UI.append('y')
            sage: dev.commit(message='modified tracked')
            Commit your changes to branch "branch1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are
            #  done.
            self._UI.info(['',
                           '(use "{0}" to checkout a branch)'
                           .format(self._format_command("checkout"))])
            self._UI.debug('Committing pending changes to branch "{0}".'.format(branch))
                    self._UI.show(['The following files in your working directory are not tracked by git:', ''] +
                                  ['    ' + f for f in untracked_files ] +
                                  [''])
                    if self._UI.confirm('Start tracking any of these files?', default=False):
                            if self._UI.confirm('Start tracking "{0}"?'.format(file), default=False):
                    from sage.dev.misc import tmp_filename
                    commit_message = tmp_filename()
                    with open(commit_message, 'w') as f:
                        f.write(COMMIT_GUIDE)
                    self._UI.edit(commit_message)
                    message = "\n".join([line for line in open(commit_message).read().splitlines() 
                                         if not line.startswith("#")]).strip()
                if not self._UI.confirm('Commit your changes to branch "{0}"?'.format(branch), default=True):
                    self._UI.info(['', 'Run "{0}" first if you want to commit to a different branch or ticket.'],
                                  self._format_command("checkout"))
                    raise OperationCancelledError("user does not want to create a commit")
                self._UI.debug("A commit has been created.")
                self._UI.info(['', 'Use "{0}" to push your commits to the trac server once you are done.'],
                              self._format_command("push"))
                self._UI.debug("Not creating a commit.")
            - :meth:`push` -- To push changes after setting the remote
              branch
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
                self._UI.error('You must specify "branch" in detached HEAD state.')
                self._UI.info(['', 'Use "{0}" to checkout a branch'],
                              self._format_command('checkout'))
                self._UI.error('no local branch for ticket #{0} found. Cannot set remote branch'
                               ' for that ticket.', ticket)
            self._UI.warning('The remote branch "{0}" is not in your user scope. You probably'
                             ' do not have permission to push to that branch.', remote_branch)
            self._UI.info(['', 'You can always use "u/{1}/{0}" as the remote branch name.'],
                          remote_branch, self.trac._username)
    def push(self, ticket=None, remote_branch=None, force=False):
        Push the current branch to the Sage repository.
          set to ``remote_branch`` after the current branch has been pushed there.
          branch to push to; if ``None``, then a default is chosen
        - ``force`` -- a boolean (default: ``False``), whether to push if
            - :meth:`commit` -- Save changes to the local repository.
            - :meth:`pull` -- Update a ticket with changes from the remote
              repository.
        TESTS:
        Alice tries to push to ticket 1 which does not exist yet::
            sage: alice.push(ticket=1)
            Ticket name "1" is not valid or ticket does not exist on trac.
        Alice creates ticket 1 and pushes some changes to it::
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
        Now Bob can check that ticket out and push changes himself::
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
        Now Alice can pull these changes::
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
        After Alice pushed her changes, Bob can not set the branch field anymore::
            sage: alice.push()
            Local commits that are not on the remote branch "u/alice/ticket/1":
            <BLANKLINE>
                ...: alice: modified tracked
                ...: bob: modified tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/ticket/1" to "u/alice/ticket/1"
            Change the "Branch:" field? [Yes/no] y
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ....: bob: added tracked2
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Not setting the branch field for ticket #1 to "u/bob/ticket/1" because
            "u/bob/ticket/1" and the current value of the branch field "u/alice/ticket/1"
            have diverged.
            <BLANKLINE>
            #  Use "sage --dev push --force --ticket=1 --remote-branch=u/bob/ticket/1" to
            #  overwrite the branch field.
            <BLANKLINE>
            #  Use "sage --dev download --ticket=1" to merge the changes introduced by the
            #  remote "u/alice/ticket/1" into your local branch.
            sage: bob.pull()
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: Merge branch 'u/alice/ticket/1' of ... into ticket/1
                ...: alice: modified tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
            sage: bob.push(2)
            Ticket name "2" is not valid or ticket does not exist on trac.
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            sage: bob.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.push(2)
            About to push the branch "ticket/1" to "u/bob/ticket/2" for ticket #2. However,
            your local branch for ticket #2 seems to be "ticket/2".
             Do you really want to proceed? [yes/No] y
            <BLANKLINE>
            #  Use "sage --dev checkout --ticket=2 --branch=ticket/1" to permanently set
            #  "ticket/1" as the branch associated to ticket #2.
            The branch "u/bob/ticket/2" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: bob.push(remote_branch="u/bob/branch1")
            The branch "u/bob/branch1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/ticket/1" to "u/bob/branch1"
            Change the "Branch:" field? [Yes/no] y
            Merging the remote branch "u/bob/ticket/2" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            sage: with open("another_file", "w") as f: f.write("bob after merge(2)")
            sage: bob.push()
            Remote branch "u/bob/ticket/1" is idential to your local branch "ticket/1
            <BLANKLINE>
            #  (use "sage --dev commit" to commit changes before pushing)
            sage: bob._UI.extend(['y', 'y', 'y'])
            sage: bob.commit(message="Bob's merge")  # oops
            The following files in your working directory are not tracked by git:
            <BLANKLINE>
                another_file
            <BLANKLINE>
            Start tracking any of these files? [yes/No] y
            Start tracking "another_file"? [yes/No] y
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are
            #  done.
            sage: bob._UI.extend(['y', 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: Bob's merge
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/branch1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
            Uploading your dependencies for ticket #1: "" => "#2"
            sage: with open("another_file", "w") as f: f.write("bob after push")
            sage: bob._UI.extend(['y', 'y', 'y'])
            sage: bob.commit(message='another commit')
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are
            #  done.
            sage: bob._UI.extend(['y', "keep", 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: another commit
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Trac ticket #1 depends on #2 while your local branch depends on no tickets.
            Updating dependencies is recommended but optional.
            Action for dependencies? [upload/download/keep] keep
            sage: with open("another_file", "w") as f: f.write("bob after 2nd push")
            sage: bob._UI.append('y')
            sage: bob.commit(message='final commit')
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are
            #  done.

            sage: bob._UI.extend(['y', 'download', 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: final commit
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Trac ticket #1 depends on #2 while your local branch depends on no tickets.
            Updating dependencies is recommended but optional.
            Action for dependencies? [upload/download/keep] download
            self._UI.error("Cannot push while in detached HEAD state.")
            raise OperationCancelledError("cannot push while in detached HEAD state")
                    raise SageDevValueError("remote_branch must be specified since #{0}"
                                            " has no remote branch set.".format(ticket))
                    raise SageDevValueError("remote_branch must be specified since the"
                                            " current branch has no remote branch set.")
        user_confirmation = False
                self._UI.show('About to push the branch "{0}" to "{1}" for ticket #{2}.'
                              ' However, your local branch for ticket #{2} seems to be "{3}".',
                              branch, remote_branch, ticket, self._local_branch_for_ticket(ticket))
                user_confirmation = self._UI.confirm(' Do you really want to proceed?', default=False)
                if user_confirmation:
                    self._UI.info(['', 
                                   'Use "{2}" to permanently set "{1}" as the branch'
                                   ' associated to ticket #{0}.'],
                                  ticket, branch, self._format_command("checkout",ticket=ticket,branch=branch))
                self._UI.show('About to push the local branch "{0}" to remote branch "{1}" for'
                              ' ticket #{2}. However, that branch is already associated to ticket #{3}.',
                              branch, remote_branch, ticket, self._ticket_for_local_branch(branch))
                user_confirmation = self._UI.confirm(' Do you really want to proceed?', default=False)
                if user_confirmation:
                    self._UI.info(['', 'Use "{2}" to permanently set the branch associated to'
                                   ' ticket #{0} to "{1}". To create a new branch from "{1}" for'
                                   ' #{0}, use "{3}" and "{4}".'], 
                                  ticket, branch, 
                                  self._format_command("checkout",ticket=ticket,branch=branch), 
                                  self._format_command("checkout",ticket=ticket), 
                                  self._format_command("merge", branch=branch))

        self._UI.debug('Pushing your changes in "{0}" to "{1}".'.format(branch, remote_branch))
                self._UI.show('The branch "{0}" does not exist on the remote server.', remote_branch)
                if not self._UI.confirm('Create new remote branch?', default=True):
                    self._UI.error('Not pushing your changes because they would discard some of'
                                   ' the commits on the remote branch "{0}".', remote_branch)
                    self._UI.info(['', 'Use "{0}" if you really want to overwrite the remote branch.'],
                                  self._format_command("push", ticket=ticket, 
                                                       remote_branch=remote_branch, force=True))
            if remote_branch_exists and not force and \
               self.git.commit_for_branch(branch) == self.git.commit_for_ref('FETCH_HEAD'):
                self._UI.show('Remote branch "{0}" is idential to your local branch "{1}',
                              remote_branch, branch)
                self._UI.info(['', '(use "{0}" to commit changes before pushing)'],
                               self._format_command("commit"))
                return
                            self._UI.show(['Local commits that are not on the remote branch "{0}":', ''] +
                                          ['    ' + c for c in commits.splitlines()] +
                                          [''], remote_branch)
                            if not self._UI.confirm('Push to remote branch?', default=True):
                    self._upload_ssh_key() # make sure that we have access to the repository
                    self.git.super_silent.push(self.git._repository, 
                                               "{0}:{1}".format(branch, remote_branch), 
                                               force=force)
            self._UI.debug('Changes in "{0}" have been pushed to "{1}".'.format(branch, remote_branch))
            self._UI.debug("Did not push any changes.")
                self._UI.debug('Not setting the branch field for ticket #{0} because it already'
                               ' points to your branch "{1}".'.format(ticket, remote_branch))
                self._UI.debug('Setting the branch field of ticket #{0} to "{1}".'.format(ticket, remote_branch))
                        self._UI.error('Not setting the branch field for ticket #{0} to "{1}" because'
                                       ' "{1}" and the current value of the branch field "{2}" have diverged.'
                                       .format(ticket, remote_branch, current_remote_branch))
                        self._UI.info(['', 'Use "{0}" to overwrite the branch field.', '',
                                       'Use "{1}" to merge the changes introduced by'
                                       ' the remote "{2}" into your local branch.'],
                                      self._format_command("push", ticket=ticket,
                                                           remote_branch=remote_branch, force=True),
                                      self._format_command("download", ticket=ticket),
                                      current_remote_branch)
                    self._UI.show('The branch field of ticket #{0} needs to be'
                                  ' updated from its current value "{1}" to "{2}"'
                                  ,ticket, current_remote_branch, remote_branch)
                    if not self._UI.confirm('Change the "Branch:" field?', default=True):
        if ticket and self._has_ticket_for_local_branch(branch):
            new_dependencies_ = self._dependencies_for_ticket(self._ticket_for_local_branch(branch))
                    self._UI.show('Trac ticket #{0} depends on {1} while your local branch depends'
                                  ' on {2}. Updating dependencies is recommended but optional.', 
                                  ticket, old_dependencies, new_dependencies or "no tickets"),
                    sel = self._UI.select('Action for dependencies?', options=("upload", "download", "keep"))
                        self._UI.debug("Setting dependencies for #{0} to {1}.", ticket, old_dependencies)
                self._UI.debug("Not uploading your dependencies for ticket #{0} because the"
                               " dependencies on trac are already up-to-date.", ticket)
                self._UI.show('Uploading your dependencies for ticket #{0}: "{1}" => "{2}"',
                              ticket, old_dependencies, new_dependencies)
    def reset_to_clean_state(self, error_unless_clean=True, helpful=True):
        - ``error_unless_clean`` -- a boolean (default: ``True``),
          whether to raise an
          :class:`user_interface_error.OperationCancelledError` if the
            sage: dev._wrap("reset_to_clean_state")
            sage: UI.append("cancel")
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] cancel
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)
            sage: UI.append("reset")
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] reset
        self._UI.show('Repository is in an unclean state ({0}).'
                      ' Resetting the state will discard any uncommited changes.',
                      ', '.join(states))
        sel = self._UI.select('Reset repository?',
                              options=('reset', 'cancel'), default=1)
        if sel == 'cancel':
            if not error_unless_clean:
            if helpful:
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
        elif sel == 'reset':
            self.git.reset_to_clean_state()
        else:
            assert False
    def clean(self, error_unless_clean=True):
        Restore the working directory to the most recent commit.
        - ``error_unless_clean`` -- a boolean (default: ``True``),
          whether to raise an
          :class:`user_interface_error.OperationCancelledError` if the
            sage: dev.clean()
        Check that nothing happens if there are only untracked files::
            sage: open("untracked","w").close()
            sage: dev.clean()
        Uncommitted changes can simply be dropped::
            sage: open("tracked","w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("discard")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            sage: dev.clean()
        Uncommitted changes can be kept::
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("cancel")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
        Or stashed::
            sage: UI.append("stash")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] stash
            Your changes have been moved to the git stash stack. To re-apply your changes
            later use "git stash apply".
            sage: dev.clean()
        """
            self.reset_to_clean_state(error_unless_clean)
        except OperationCancelledError:
            self._UI.error("Can not clean the working directory unless in a clean state.")
            raise

        if not self.git.has_uncommitted_changes():
            return
        files = [line[2:] for line in self.git.status(porcelain=True).splitlines()
                 if not line.startswith('?')]

        self._UI.show(
            ['The following files in your working directory contain uncommitted changes:'] +
            [''] +
            ['    ' + f for f in files ] +
            [''])
        cancel = 'cancel' if error_unless_clean else 'keep'
        sel = self._UI.select('Discard changes?',
                              options=('discard', cancel, 'stash'), default=1)
        if sel == 'discard':
            self.git.clean_wrapper()
        elif sel == cancel:
            if error_unless_clean:
                raise OperationCancelledError("User requested not to clean the working directory.")
        elif sel == 'stash':
            self.git.super_silent.stash()
            self._UI.show('Your changes have been moved to the git stash stack. '
                          'To re-apply your changes later use "git stash apply".')
        else:
            assert False
            :meth:`create_ticket`, :meth:`comment`,
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
    def needs_review(self, ticket=None, comment=''):
        r"""
        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.
            :meth:`set_positive_review`, :meth:`comment`,
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.needs_review(comment='Review my ticket!')
        self._UI.debug("Ticket #%s marked as needing review"%ticket)
    def needs_work(self, ticket=None, comment=''):
        r"""
        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.
            :meth:`set_positive_review`, :meth:`comment`,
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.needs_work(comment='Need to add an untracked file!')
        self._UI.debug("Ticket #%s marked as needing work"%ticket)
    def needs_info(self, ticket=None, comment=''):
        r"""
        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.
            :meth:`edit_ticket`, :meth:`needs_review`,
            :meth:`positive_review`, :meth:`comment`,
            :meth:`needs_work`
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.needs_info(comment='Why is a tracked file enough?')
        self._UI.debug("Ticket #%s marked as needing info"%ticket)
    def positive_review(self, ticket=None, comment=''):
        r"""
        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.
            :meth:`edit_ticket`, :meth:`needs_review`,
            :meth:`needs_info`, :meth:`comment`,
            :meth:`needs_work`
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.positive_review()
        self._UI.debug("Ticket #%s reviewed!"%ticket)
    def comment(self, ticket=None):
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.comment()
            :meth:`edit_ticket`, :meth:`comment`,
            ticket must be specified if not currently on a ticket.
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            Your branch "ticket/1" has 0 commits.
        After pushing the local branch::
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Your branch "ticket/1" has 0 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 0 commits. It does not differ from "ticket/1".
            Your branch "ticket/1" has 1 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 0 commits. "ticket/1" is ahead of "u/doctest/ticket/1" by 1 commits:
        Pushing them::
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/1":
            <BLANKLINE>
                ...: added tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Your branch "ticket/1" has 1 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 1 commits. It does not differ from "ticket/1".
            Your branch "ticket/1" has 0 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 1 commits. "u/doctest/ticket/1" is ahead of "ticket/1" by 1 commits:
            sage: dev.push(remote_branch="u/doctest/branch1", force=True)
            The branch "u/doctest/branch1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Your branch "ticket/1" has 2 commits.
            The trac ticket points to the branch "u/doctest/branch1" which has 3 commits. "u/doctest/branch1" is ahead of "ticket/1" by 1 commits:
            Your remote branch "u/doctest/ticket/1" has 1 commits. The branches "u/doctest/ticket/1" and "ticket/1" have diverged.
            "u/doctest/ticket/1" is ahead of "ticket/1" by 1 commits:
            "ticket/1" is ahead of "u/doctest/ticket/1" by 2 commits:
        commits = lambda a, b: list(reversed(
            self.git.log("{0}..{1}".format(a,b), "--pretty=%an <%ae>: %s").splitlines()))
                return 'It does not differ from "{0}".'.format(b)
                return '"{0}" is ahead of "{1}" by {2} commits:\n{3}'.format(a,b,len(b_to_a), "\n".join(b_to_a))
                return '"{0}" is ahead of "{1}" by {2} commits:\n{3}'.format(b,a,len(a_to_b),"\n".join(a_to_b))
                return ('The branches "{0}" and "{1}" have diverged.\n"{0}" is ahead of'
                        ' "{1}" by {2} commits:\n{3}\n"{1}" is ahead of "{0}" by {4}'
                        ' commits:\n{5}'.format(a, b, len(b_to_a), "\n".join(b_to_a), 
                                                len(a_to_b), "\n".join(a_to_b)))
            local_summary = 'Your branch "{0}" has {1} commits.'.format(branch, len(master_to_branch))
                ticket_summary = 'The trac ticket points to the branch "{0}" which does not exist.'
                ticket_summary = 'The trac ticket points to the' \
                    ' branch "{0}" which has {1} commits.'.format(ticket_branch, len(master_to_ticket))
                        ticket_summary += ' The branch can not be compared to your local' \
                            ' branch "{0}" because the branches are based on different versions' \
                            ' of sage (i.e. the "master" branch).'
            remote_summary = 'Your remote branch "{0}" has {1} commits.'.format(
                remote_branch, len(master_to_remote))
                    remote_summary += ' The branch can not be compared to your local' \
                        ' branch "{0}" because the branches are based on different version' \
                        ' of sage (i.e. the "master" branch).'
    def prune_tickets(self):
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.tickets()
                : master
            * #1: ticket/1 summary
            sage: dev.prune_tickets()
            sage: dev.tickets()
                : master
            * #1: ticket/1 summary
            sage: dev.prune_tickets()
            Cannot delete "ticket/1": is the current branch.
            <BLANKLINE>
            #  (use "sage --dev vanilla" to switch to the master branch)
            sage: dev.prune_tickets()
            Moved your branch "ticket/1" to "trash/ticket/1".
            sage: dev.tickets()
            sage: dev.prune_tickets()
                    self.abandon(ticket, helpful=False)
    def abandon(self, ticket_or_branch=None, helpful=True):
        - ``helpful`` -- boolean (default: ``True``). Whether to print
          informational messages to guide new users.

            - :meth:`prune_tickets` -- abandon tickets that have
              been closed.
            - :meth:`tickets` -- list local non-abandoned tickets.
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Cannot delete "ticket/1": is the current branch.
            <BLANKLINE>
            #  (use "sage --dev vanilla" to switch to the master branch)
            Moved your branch "ticket/1" to "trash/ticket/1".
            <BLANKLINE>
            #  Use "sage --dev checkout --ticket=1 --base=master" to restart working on #1
            #  with a clean copy of the master branch.
            sage: dev.checkout(ticket=1, base=MASTER_BRANCH)
            About to create a new branch for #1 based on "master". However, the trac ticket
            for #1 already refers to the branch "u/doctest/ticket/1". The new branch will
            not contain any work that has already been done on "u/doctest/ticket/1".
            Create fresh branch? [yes/No] y
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
                raise SageDevValueError("Cannot abandon #{0}: no local branch for this ticket.", ticket)
                self._UI.error("Cannot delete the master branch.")
                    self._UI.error('Cannot delete "{0}": is the current branch.', branch)
                    self._UI.info(['', '(use "{0}" to switch to the master branch)'],
                                  self._format_command("vanilla"))
            self._UI.show('Moved your branch "{0}" to "{1}".', branch, new_branch)
            if helpful:
                self._UI.info(['',
                               'Use "{0}" to restart working on #{1} with a clean copy of the master branch.'], 
                               self._format_command("checkout", ticket=ticket, base=MASTER_BRANCH), ticket)
            - :meth:`merge` -- merge into the current branch rather
              than creating a new one
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            self.clean()
        self._UI.debug('Creating a new branch "{0}".'.format(branch))
                self._UI.debug('Merging {2} branch "{0}" into "{1}".'
                              .format(branch_name, branch, local_remote))
                self.merge(branch, pull=local_remote=="remote")
            self.git.clean_wrapper()
            self._UI.debug('Deleted branch "{0}".'.format(branch))
    def merge(self, ticket_or_branch=MASTER_BRANCH, pull=None, create_dependency=None):
          ticket, if ``pull`` is ``False``), for the name of a local or
        - ``pull`` -- a boolean or ``None`` (default: ``None``); if
          ``ticket_or_branch`` identifies a ticket, whether to pull the
          ``ticket_or_branch`` is a branch name, then ``pull`` controls
            the remote server during :meth:`push` and :meth:`pull`.
            - :meth:`show_dependencies` -- see the current
              dependencies.
            - :meth:`GitInterface.merge` -- git's merge command has
              more options and can merge multiple branches at once.
            - :meth:`gather` -- creates a new branch to merge into
              rather than merging into the current branch.
        TESTS:
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            Cannot merge remote branch for #1 because no branch has been set on the trac
            ticket.
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: alice.merge("#1", pull=False)
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
        A remote branch for a local branch is only merged in if ``pull`` is set::
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: alice.merge("ticket/1", pull=True)
            Branch "ticket/1" does not exist on the remote system.
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/2".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort
            sage: alice._UI.append("ok")
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/2".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] ok
            Created a commit from your conflict resolution.
            cannot merge a ticket into itself
        We also cannot merge if the working directory has uncommited changes::

            sage: alice._UI.append("cancel")
            sage: with open("alice2","w") as f: f.write("uncommited change")
            sage: alice.merge(1)
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 alice2
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Cannot merge because working directory is not in a clean state.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your changes)
            self.clean()
            self._UI.info(['', '(use "{0}" to commit your changes)'], 
                          self._format_command('commit'))
            self._UI.error('Not on any branch.')
            self._UI.info(['', '(use "{0}" to checkout a branch)'],
                           self._format_command("checkout"))
            if pull == False:
                raise SageDevValueError('"pull" must not be "False" when merging dependencies.')
                raise SageDevValueError('"create_dependency" must not be set when merging dependencies.')
                self._UI.debug("Merging dependency #{0}.".format(dependency))
                self.merge(ticket_or_branch=dependency, pull=True)
            if pull is None:
                pull = True
            if pull:
                    self._UI.error("Cannot merge remote branch for #{0} because no branch has"
                                   " been set on the trac ticket.", ticket)
        elif pull == False or (pull is None and not 
                               self._is_remote_branch_name(ticket_or_branch, exists=True)):
            pull = False
                    raise SageDevValueError('"create_dependency" must not be "True" if'
                                            ' "ticket_or_branch" is a local branch which'
                                            ' is not associated to a ticket.')
            pull = True
                raise SageDevValueError('"create_dependency" must not be "True" if'
                                        ' "ticket_or_branch" is a local branch.')
        if pull:
                self._UI.error('Can not merge remote branch "{0}". It does not exist.',
                               remote_branch)
            self._UI.show('Merging the remote branch "{0}" into the local branch "{1}".',
                          remote_branch, current_branch)
            self._UI.show('Merging the local branch "{0}" into the local branch "{1}".',
                          branch, current_branch)
            self._UI.show('Automatic merge successful.')
            self._UI.info(['', '(use "{0}" to commit your merge)'], 
                          self._format_command('commit'))
                self._UI.show('Automatic merge failed, there are conflicting commits.')
                excluded = ['Aborting', 
                    "Automatic merge failed; fix conflicts and then commit the result."]
                lines = [line for line in lines if line not in excluded]
                self._UI.show([''] + lines + [''])
                self._UI.show('Please edit the affected files to resolve the conflicts.'
                              ' When you are finished, your resolution will be commited.')
                sel = self._UI.select("Finished?", ['ok', 'abort'], default=1)
                if sel == 'ok':
                    self._UI.show("Created a commit from your conflict resolution.")
                elif sel == 'abort':
                else:
                    assert False
                self.git.clean_wrapper()
                self._UI.debug("Not recording dependency on #{0} because #{1} already depends on #{0}.",
                               ticket, current_ticket)
                self._UI.show(['', "Added dependency on #{0} to #{1}."], ticket, current_ticket)
    def tickets(self, include_abandoned=False, cached=True):
        - ``cached`` -- boolean (default: ``True``), whether to try to pull the
          summaries from the ticket cache; if ``True``, then the summaries
          might not be accurate if they changed since they were last updated.
          To update the summaries, set this to ``False``.

            - :meth:`abandon_ticket` -- hide tickets from this method.
            - :meth:`remote_status` -- also show status compared to
              the trac server.
            - :meth:`current_ticket` -- get the current ticket.
            sage: dev.tickets()
            * : master
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.tickets()
                : master
              #1: ticket/1 summary
            * #2: ticket/2 summary
        from git_error import DetachedHeadError
        try:
            current_branch = self.git.current_branch()
        except DetachedHeadError:
            current_branch = None
        ret = []
        for branch in branches:
            ticket = None
            ticket_summary = ""
            extra = " "
            if self._has_ticket_for_local_branch(branch):
                ticket = self._ticket_for_local_branch(branch)
                try:
                    try:
                        ticket_summary = self.trac._get_attributes(ticket, cached=cached)['summary']
                    except KeyError:
                        ticket_summary = self.trac._get_attributes(ticket, cached=False)['summary']
                except TracConnectionError:
                    ticket_summary = ""
            if current_branch == branch:
                extra = "*"
            ticket_str = "#"+str(ticket) if ticket else ""
            ret.append(("{0:>7}: {1} {2}".format(ticket_str, branch, ticket_summary), extra))
        while all([info.startswith(' ') for (info, extra) in ret]):
            ret = [(info[1:],extra) for (info, extra) in ret]
        ret = sorted(ret)
        ret = ["{0} {1}".format(extra,info) for (info,extra) in ret]
        self._UI.show("\n".join(ret))

    def vanilla(self, release=MASTER_BRANCH):
        Return to a clean version of Sage.
        - ``release`` -- a string or decimal giving the release name (default:
          ``'master'``).  In fact, any tag, commit or branch will work.  If the
          tag does not exist locally an attempt to fetch it from the server
          will be made.
            - :meth:`checkout` -- checkout another branch, ready to
              develop on it.
            - :meth:`pull` -- pull a branch from the server and merge
              it.
        release = str(release)
            self.clean()
            self._UI.error("Cannot checkout a release while your working directory is not clean.")
                self._UI.error('"{0}" does not exist locally or on the remote server.'.format(release))
            - :meth:`commit` -- record changes into the repository.
            - :meth:`tickets` -- list local tickets (you may
              want to commit your changes to a branch other than the
              current one).
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/2" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Created ticket #3 at https://trac.sagemath.org/3.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=3" to create a new local branch)
            sage: dev.checkout(ticket=3)
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            The branch "u/doctest/ticket/3" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            Merging the remote branch "u/doctest/ticket/1" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Merging the remote branch "u/doctest/ticket/2" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            sage: dev.checkout(ticket="#1")
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.checkout(ticket="#2")
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/2":
            <BLANKLINE>
                ...: added ticket2
            <BLANKLINE>
            Push to remote branch? [Yes/no] y

            sage: dev.checkout(ticket="#3")
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev commit" to save changes in a new commit when you are finished
            #  editing.
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/3":
            <BLANKLINE>
                ...: added ticket3
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Uploading your dependencies for ticket #3: "" => "#1, #2"