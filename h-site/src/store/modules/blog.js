import db from '../../firebaseConfig.js'

const state = {
	posts: [],
	isFetched: false
}

const getters = {
	allPosts: (state) => state.posts
}

const actions = {
	async getPosts({ commit }) {
		if (!state.isFetched) {
			await db.blog.get().then(
			function(querySnapshot) {
				querySnapshot.forEach(function(doc) {
					//console.log(doc.id, "=>", doc.data())
					commit('newPost', {
						id: doc.id,
						title: doc.data().title,
						post: doc.data().post,
						likes: doc.data().likes,
						timeOfPost: doc.data().timeOfPost
					})
				})
				state.isFetched = !state.isFetched
				//commit('setTodos', querySnapshot)
				//console.log("Todos retrieved from database: ", state.isFetched)
			})
		}
	},
	async newLike(updpost) {
		console.log(updpost)
		await db.blog.doc(updpost).update({
			likes: db.incerement
		}).then(function() {
			console.log("like added")
		}).catch(function(error) {
			console.error("An error occured liking a post: ", error)
		})
	},
	async commitPost({ commit }, postToCommit) {
		await db.blog.add({
			title: postToCommit.title,
			post: postToCommit.postText,
			likes: 0,
			timeOfPost: new Date()
		}).then(function(doc) {
			console.log("Post added")
			commit('newPost', {
				id: doc.id,
				title: doc.data().title,
				post: doc.data().post,
				likes: doc.data().likes,
				timeOfPost: doc.data().timeOfPost
			})
		}).catch(function(error) {
			console.error("An error occured adding the new post: ", error)
		})
	}
}

const mutations = {
	newPost: (state, post) => state.posts.unshift(post)
}

export default {
	state,
	getters,
	actions,
	mutations
}